import os
import re
import sys
from collections import defaultdict
from glob import glob
import numpy as np
from scipy import io as sio
from scipy import sparse as sp


# Quantities stored as integer cell indices (labels).
INT_QUANTS = {'lowerAddr', 'upperAddr', 'interfaceRow', 'interfaceCol'}
# Quantities read per solve. The matrix values (diag/lower/interfaceLower) and the
# addressing (lowerAddr/upperAddr/interfaceRow/interfaceCol) are deduplicated on the
# .H side against the most recent dump: a file is written only when it changes,
# so an absent file means "unchanged since the previous dump" and the most
# recently loaded value is carried forward. (b/sol/init are written every solve.)
LOAD_QUANTS = ['diag', 'lower', 'lowerAddr', 'upperAddr',
               'interfaceRow', 'interfaceCol', 'interfaceLower', 'b']


def _parse_block(buf, dtype):
    """Parse a binary OpenFOAM `<count>\\n(<data>)` payload. The data begins at
    the first '(' (neither a FoamFile header nor the count contains one) and the
    count is the integer immediately before it. An empty field is written as a
    bare '0' with no parens, which is the p < 0 case."""
    p = buf.find(b'(')
    if p < 0:
        return np.empty(0, dtype=dtype)
    count = int(re.search(rb'(\d+)\s*$', buf[:p]).group(1))
    return np.frombuffer(buf, dtype=dtype, count=count, offset=p + 1).copy()


def _load_uncollated(fname, dtype):
    """Read a plain (single-processor) OpenFOAM IOField/IOList."""
    return _parse_block(open(fname, 'rb').read(), dtype)


def _load_collated(fname, dtype, block):
    """Extract one processor block from an OpenFOAM collated (decomposedBlockData)
    file. Blocks are located by their '// processorN <nbytes> (' markers, using
    the byte-length prefixes to skip the (binary) payloads. The master header has
    no such marker and the '// * * *' separator is not always present, so the scan
    starts from the beginning of the file."""
    data = open(fname, 'rb').read()
    pat = re.compile(rb'//\s*processor\d+\s+(\d+)\s*\(')
    pos = 0
    spans = []
    while True:
        mo = pat.search(data, pos)
        if not mo:
            break
        nbytes = int(mo.group(1))
        start = mo.end()
        spans.append((start, nbytes))
        pos = start + nbytes
    start, nbytes = spans[block]
    return _parse_block(data[start:start + nbytes], dtype)


def load(fname, dtype=np.float64, block=None):
    if block is None:
        return _load_uncollated(fname, dtype)
    return _load_collated(fname, dtype, block)


def read_info(fname):
    """Read the nIterations entry from a solve's info dict (the collated info is
    a single global dict, so no per-block handling is required)."""
    info = {}
    with open(fname, errors='replace') as f:
        for raw in f:
            if not raw[0] in {'/', ' '}:
                tok = raw.split()
                if len(tok) == 2:
                    entry = float(tok[1].rstrip(';'))
                    if not tok[0][-4:] in {'ance', 'dual', 'Time'}:
                        entry = int(entry)
                    info[tok[0]] = entry
    return info


def openfoam_residual(A, b, x):

    Ax = A @ x
    Axbar = A @ np.ones(x.shape) * x.mean()
    return np.absolute(b - Ax).sum() / (np.absolute(Ax-Axbar) + np.absolute(b-Axbar)).sum()


def discover(root):
    """Map (field, timestep, solve) -> solve directory for one subdomain root.

    Layout: <root>/<timestep>/<field>-Absol/<solve>/<files>, with an extra
    nesting level (e.g. <timestep>/<region>/<field>-Absol/...) also handled.
    """
    out = {}
    for tdir in glob(f'{root}/*'):
        try:
            timestep = float(tdir.rsplit('/', 1)[-1])
        except ValueError:
            continue
        for depth in (1, 2):
            for fdir in glob(tdir + '/*' * depth):
                if not os.path.isdir(fdir):
                    continue
                field = '.'.join(fdir.split('/')[-depth:])
                for sdir in glob(f'{fdir}/*'):
                    try:
                        solve = int(sdir.rsplit('/', 1)[-1])
                    except ValueError:
                        continue
                    out[(field, timestep, solve)] = sdir
    return out


def subdomains(folder):
    """Return [(root, block), ...]. block is an int for a collated processorsN
    directory (one entry per processor block), else None (uncollated processorN
    directories, or a serial / reconstructed tree)."""
    coll = sorted(d for d in glob(f'{folder}/processors*')
                  if re.fullmatch(r'processors\d+', os.path.basename(d)))
    if coll:
        specs = []
        for root in coll:
            n = int(re.fullmatch(r'processors(\d+)', os.path.basename(root)).group(1))
            specs += [(root, b) for b in range(n)]
        return specs
    unc = sorted((d for d in glob(f'{folder}/processor*')
                  if re.fullmatch(r'processor\d+', os.path.basename(d))),
                 key=lambda d: int(os.path.basename(d)[len('processor'):]))
    if unc:
        return [(d, None) for d in unc]
    return [(folder, None)]


def build_orig(specs, region, sizes, offs, nGlobal):
    """Permutation from the decomposed (globalIndex) cell order to the original,
    undecomposed cell order, read from each subdomain's cellProcAddressing. So the
    assembled matrix carries the original mesh numbering and a parallel run matches
    a serial one element-wise. Returns None when the addressing is absent (serial /
    reconstructed tree), i.e. the ordering is already the original one."""
    orig = np.empty(nGlobal, dtype=np.int64)
    for si, (root, block) in enumerate(specs):
        cpa = os.path.join(root, 'constant', *region, 'polyMesh', 'cellProcAddressing')
        if not os.path.isfile(cpa):
            return None
        a = load(cpa, np.uint32, block).astype(np.int64)   # proc-local cell -> original cell
        if a.shape[0] != sizes[si]:
            return None
        orig[offs[si]:offs[si + 1]] = a
    return orig


if __name__ == '__main__':

    folder = sys.argv[1].rstrip('/')
    specs = subdomains(folder)
    maps = {}
    for root, _ in specs:
        maps.setdefault(root, discover(root))
    allkeys = sorted(set().union(*(set(maps[r]) for r, _ in specs)))

    state = defaultdict(dict)    # subdomain index -> {quantity: array}
    counter = defaultdict(int)   # field -> .mat index
    orig_cache = {}              # region -> permutation (globalIndex order -> original cell)

    mrr = 0.0
    ml2 = 0.0

    for n, key in enumerate(allkeys):
        field, timestep, solve = key

        # --- read every subdomain block, carrying forward absent (addressing) data
        blocks = []
        sds = []
        complete = True
        for si, (root, block) in enumerate(specs):
            sd = maps[root].get(key)
            if sd is None or not os.path.isfile(f'{sd}/info'):
                complete = False
                break
            st = state[si]
            for q in LOAD_QUANTS:
                path = f'{sd}/{q}'
                if os.path.isfile(path):
                    st[q] = load(path, np.uint32 if q in INT_QUANTS else np.float64, block)
            info = read_info(f'{sd}/info')
            sol_b = load(f'{sd}/sol', np.float64, block)
            try:
                init_b = load(f'{sd}/init', np.float64, block)
            except (ValueError, OSError):
                init_b = None
            blocks.append((info, dict(st), sol_b, init_b))
            sds.append(sd)
        if not complete:
            continue
        if any('lowerAddr' not in st for _, st, _, _ in blocks):
            print(f"\nskip {key}: addressing not seen yet (retain the run's first dump)")
            continue

        # --- per-subdomain offsets in the decomposed (globalIndex) order, which
        #     is how the dumped global addressing is numbered
        sizes = [st['diag'].shape[0] for _, st, _, _ in blocks]
        offs = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
        nGlobal = int(offs[-1])
        niter = blocks[0][0].get('nIterations', 0)

        # --- map that order back to the original undecomposed cell numbering via
        #     cellProcAddressing (cached per region; identity when not decomposed),
        #     so the output matches a serial run element-wise
        region = tuple(os.path.relpath(sds[0], specs[0][0]).split('/')[1:-2])
        if region not in orig_cache:
            orig_cache[region] = build_orig(specs, region, sizes, offs, nGlobal)
        orig = orig_cache[region]
        if orig is None:
            orig = np.arange(nGlobal, dtype=np.int64)

        # --- assemble the monolithic global symmetric matrix M (original ordering)
        rows, cols, vals = [], [], []
        bvec = np.zeros(nGlobal)
        solvec = np.zeros(nGlobal)
        initvec = np.zeros(nGlobal)
        for si, (info, st, sol_b, init_b) in enumerate(blocks):
            go = orig[offs[si]:offs[si + 1]]                             # this block's original cell ids
            rows += [go, orig[st['upperAddr']], orig[st['lowerAddr']]]   # diagonal + internal
            cols += [go, orig[st['lowerAddr']], orig[st['upperAddr']]]
            vals += [st['diag'], st['lower'], st['lower']]
            ir = st.get('interfaceRow')                                      # interface couplings
            if ir is not None and ir.size:
                rows.append(orig[ir]); cols.append(orig[st['interfaceCol']]); vals.append(st['interfaceLower'])
            bvec[go] = st['b']
            solvec[go] = sol_b
            if init_b is not None and init_b.size:
                initvec[go] = init_b

        M = sp.coo_array(
            (np.concatenate(vals).astype(np.float64),
             (np.concatenate(rows).astype(np.int64),
              np.concatenate(cols).astype(np.int64))),
            shape=(nGlobal, nGlobal),
        ).tocsr().tocoo()   # tocsr() sums any coincident entries

        # --- de-symmetrization (storage only): keep the diagonal and the
        #     doubled strict upper triangle; the true matrix is 0.5*(A + A.T).
        dm = M.row == M.col
        um = M.row < M.col
        A = sp.coo_array(
            (np.concatenate([M.data[dm], 2.0 * M.data[um]]),
             (np.concatenate([M.row[dm], M.row[um]]),
              np.concatenate([M.col[dm], M.col[um]]))),
            shape=(nGlobal, nGlobal),
        )

        outdir = f'{folder}/{field}-systems'
        os.makedirs(outdir, exist_ok=True)
        counter[field] += 1
        sio.savemat(
            os.path.join(outdir, f'{counter[field]}.mat'),
            {'A': A, 'b': bvec, 'sol': solvec, 'init': initvec, 'niter': niter},
        )

        print('processed', n + 1, '/', len(allkeys), 'systems', end='\t')
        symA = 0.5 * (A + A.T)
        effTol = max(blocks[0][0]['relativeTolerance'] * openfoam_residual(symA, bvec, initvec),
                     blocks[0][0]['tolerance'])

        mrr = max(mrr, openfoam_residual(symA, bvec, solvec) / effTol)
        print('max residual:reported:', round(mrr, 2), end='\t')
        ml2 = max(ml2, np.linalg.norm(symA @ solvec - bvec) / np.linalg.norm(bvec))
        print('max l2 error:', ml2, end='\r')

    print()
