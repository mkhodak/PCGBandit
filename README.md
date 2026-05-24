# PCGBandit

Package for tuning PCG preconditioners in OpenFOAM simulations on-the-fly.
Includes an implementation of a thresholded incomplete Cholesky (ICT) preconditioner.
This repository is for the publication version of the software; an updated version with additional features is maintainted [here](https://github.com/the-lens-project/PCGBandit).

## Usage

Assuming only a Docker installation, an OpenCFD container with the compiled code can be launched by running `sh launch.sh`.
Assuming an existing compatible OpenFOAM installation, PCGBandit can be compiled by executing `wmake libso` in the `src/PCGBandit` and `src/ICTC` directories.
Once either is complete, to try PCGBandit go to your OpenFOAM case directory's `system` subfolder and do the following:

1. Add the line `libs ( libICTCPreconditioner.so libPCGBandit.so );` to the  `controlDict` file  
2. For any PCG solver to be tuned, replace its `solver` and `preconditioner` specifications in the `fvSolution` file with (e.g.) `solver PCGBandit; preconditioner separate; smootherTune yes; nCellsInCoarsestLevelTune yes; mergeLevelsTune yes; numDroptols 8;`. Keywords ending with `Tune` indicate tuning those keywords of the `GAMG` preconditioner. The `numDroptols` keyword indicates how many thresholds to consider when tuning the `ICTC` preconditioner.

The `ICTC` preconditioner may also be used on its own as a preconditioner for `PCG` by replacing (e.g.) `preconditioner DIC;` with `preconditioner ICTC;` in the relevant `fvSolution` file.
Similarly, PCGBandit can also be run without ICT by setting `numDroptols 0;` (the default).

## Examples

The first argument to the script specifies the simulation while the second specifies the number of cores.
Below are the commands that can be executed inside the Docker container launched by running `sh launch.sh`, or alternatively assuming the variable `$FOAM` is set to the home directory of a compatible OpenFOAM installation.
The last two assume `examples/FreeMHD.zip` has been unzipped.
```
bash run.sh boxTurb32           # OpenFOAM tutorial DNS/dnsFoam/boxTurb16 (at 2x resolution)
bash run.sh pitzDaily           # OpenFOAM tutorial incompressible/pimpleFoam/RAS/pitzDaily (at 2x resolution)
bash run.sh interStefanProblem  # OpenFOAM tutorial verificationAndValidation/multiphase/StefanProblem (at 2x resolution)
bash run.sh porousDamBreak      # OpenFOAM tutorial verificationAndValidation/multiphase/interIsoFoam/porousDamBreak (at 2x resolution)
bash run.sh closedPipe 16       # FreeMHD case file
bash run.sh fringingBField 16   # FreeMHD case file
```

## References

1. Khodak, Jung, Wynne, Chow, Kolemen. *PCGBandit: One-shot acceleration of transient PDE solvers via online-learned preconditioners.* 2025.
2. Khodak, Chow, Balcan, Talwalkar. *Learning to relax: Setting solver parameters across a sequence of linear system instances.* ICLR 2024.
3. Wynne, Saenz, Al-Salami, Xu, Sun, Hu, Hanada, Kolemen. *FreeMHD: Validation and verification of the open-source, multi-domain, multi-phase solver for electrically conductive flows.* Phys. Plasmas **32** (1) 2025.
