# PCGBandit

Package for tuning PCG preconditioners in OpenFOAM simulations on-the-fly.
Includes an implementation of a thresholded incomplete Cholesky (ICT) preconditioner.
Everything in this package assumes an existing OpenFOAM installation.

## Usage

To try PCGBandit do the following:

1. `cd src/PCGBandit && wmake libso && cd ../..`
2. `cd src/ICTC && wmake libso && cd ../..`
3. in your OpenFOAM case directory's `system` subfolder:  
    * add the line `libs ( libICTCPreconditioner.so libPCGBandit.so );` to the  `controlDict` file  
    * for any PCG solver to be tuned, replace its `solver` and `preconditioner` specifications in the `fvSolution` file with (e.g.) `solver PCGBandit; preconditioner separate; smootherTune yes; nCellsInCoarsestLevelTune yes; mergeLevelsTune yes; numDroptols 8;`. Keywords ending with `Tune` indicate tuning those keywords of the `GAMG` preconditioner. The `numDroptols` keyword indicates how many thresholds to consider when tuning the `ICTC` preconditioner.

The `ICTC` preconditioner may also be used on its own as a preconditioner for `PCG` by replacing (e.g.) `preconditioner DIC;` with `preconditioner ICTC;` in the relevant `fvSolution` file.
Similarly, PCGBandit can also be run without ICT by setting `numDroptols 0;` (the default).

## Examples

Assuming both `src/PCGBandit` and `src/ICTC` have been compiled, the script `examples/run.sh` can be used.
The first argument to the script specifies the simulation while the second specifies the number of cores.
Below are some example commands (only these cases are currently implemented in `run.sh`).
The last two assume the [FreeMHD solver](https://github.com/PlasmaControl/FreeMHD/tree/main/MHD_Solvers/solvers/epotMultiRegionInterFoam) `epotMultiRegionInterFoam` has been compiled.
```
sh run.sh boxTurb32           # OpenFOAM tutorial DNS/dnsFoam/boxTurb16 (at 2x resolution)
sh run.sh pitzDaily           # OpenFOAM tutorial incompressible/pimpleFoam/RAS/pitzDaily (at 2x resolution)
sh run.sh interStefanProblem  # OpenFOAM tutorial verificationAndValidation/multiphase/StefanProblem (at 2x resolution)
sh run.sh porousDamBreak      # OpenFOAM tutorial verificationAndValidation/multiphase/interIsoFoam/porousDamBreak (at 2x resolution)
sh run.sh closedPipe 16       # FreeMHD case file
sh run.sh fringingBField 16   # FreeMHD case file
```

## References

1. Khodak, Jung, Wynne, Chow, Kolemen. *PCGBandit: One-shot acceleration of transient PDE solvers via online-learned preconditioners.* 2025.
2. Khodak, Chow, Balcan, Talwalkar. *Learning to relax: Setting solver parameters across a sequence of linear system instances.* ICLR 2024.
3. Wynne, Saenz, Al-Salami, Xu, Sun, Hu, Hanada, Kolemen. *FreeMHD: Validation and verification of the open-source, multi-domain, multi-phase solver for electrically conductive flows.* Phys. Plasmas **32** (1) 2025.
