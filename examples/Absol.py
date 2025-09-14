import os
import sys
from collections import defaultdict
from glob import glob
from operator import itemgetter
import numpy as np
from scipy import io as sio
from scipy import sparse as sp


def getdirs(folder):
    return sorted(((int(directory.split('/')[-1]), directory) 
                   for directory in glob(f'{folder}/*')), 
                  key=itemgetter(1))


def load(fname, dtype=np.float64):
    try:
        return np.loadtxt(fname, dtype=dtype, skiprows=21, max_rows=int(np.loadtxt(fname, skiprows=19, max_rows=1)))
    except UnicodeDecodeError:
        with open(fname, 'rb') as f:
            i = 0
            while True:
                line = f.readline()
                if i or '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *' in str(line):
                    i += 1
                if i == 4:
                    break
            f.read(1)
            return np.fromfile(f, dtype=dtype)[:int(line)]


if __name__ == '__main__':

    folder = sys.argv[1]
    timesteps = []
    directories = []
    for directory in glob(f'{folder}/*'):
        try:
            timesteps.append(float(directory.split('/')[-1]))
            directories.append(directory)
        except ValueError:
            pass

    fields = defaultdict(lambda: defaultdict(lambda: 0))
    for depth in range(1, 3):
        for i, (timestep, directory) in enumerate(sorted(zip(timesteps, directories), key=itemgetter(0))):
            for subdir in glob(f'{directory}' + '/*'*depth):
                field = '.'.join(subdir.split('/')[-depth:])
                if os.path.isdir(subdir):
                    solves = []
                    solvedirs = []
                    for directory in glob(f'{subdir}/*'):
                        try:
                            solves.append(int(directory.split('/')[-1]))
                            solvedirs.append(directory)
                        except ValueError:
                            pass
                    for solve, solvedir in sorted(zip(solves, solvedirs), key=itemgetter(0)):
                        failed = False
                        for q in ['diag', 'lowerAddr', 'upperAddr', 'lower', 'A1', 'b']:
                            if os.path.isfile(f'{solvedir}/{q}'):
                                try:
                                    fields[field][q] = load(f'{solvedir}/{q}', dtype=np.uint32 if 'Addr' in q else np.float64)
                                except ValueError:
                                    failed = True
                                    break
                        if not failed:
                            b = fields[field]['b']
                            n = b.shape[0]
                            lowerAddr = fields[field]['lowerAddr']
                            upperAddr = fields[field]['upperAddr']
                            lower = fields[field]['lower']
                            diag = fields[field]['diag']
                            r = np.concatenate([np.arange(n), lowerAddr, upperAddr])
                            c = np.concatenate([np.arange(n), upperAddr, lowerAddr])
                            v = np.concatenate([diag, lower, lower])
                            A = sp.coo_array((v, (r, c)))
                            v = np.concatenate([diag + fields[field]['A1'] - A.dot(np.ones(n)), 2.0*lower])
                            A = sp.coo_array((v, (r[:v.shape[0]], c[:v.shape[0]])))
                            sol = load(f'{solvedir}/sol', dtype=np.float64)
                            try:
                                init = load(f'{solvedir}/init', dtype=np.float64)
                            except ValueError:
                                init = np.zeros(n, dtype=np.float64)
                            with open(f'{solvedir}/info', 'r') as f:
                                for line in f:
                                    if line[:11] == 'nIterations':
                                        niter = int(line.split()[1][:-1])
                                        break
                            outdir = f'{folder}/{field}-systems'
                            os.makedirs(outdir, exist_ok=True)
                            fields[field]['system'] += 1
                            sio.savemat(f'{outdir}/{fields[field]['system']}.mat', {'A': A, 'b': b, 'sol': sol, 'init': init, 'niter': niter})
                            #print(np.linalg.norm(0.5*(A+A.T).dot(sol)-b) / np.linalg.norm(b))

            print("processed", i+1, "/", len(timesteps), "timesteps", end='\r')

    print()
