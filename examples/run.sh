#!/bin/bash
#SBATCH --job-name=pitzDaily
#SBATCH --array=[0-2]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --nodelist=della-r3c[1-4]n[1-16]
#SBATCH --nodes=1

NAME=$SLURM_JOB_NAME
NPROC=$SLURM_NTASKS
SEED=$SLURM_ARRAY_TASK_ID
if [ -z $NAME ] ; then
  NAME=$1
  if [ -z $NAME ] ; then
    NAME=pitzDaily
  fi
  NPROC=$2
  if [ -z $NPROC ] ; then
    NPROC=1  
  fi
  SEED=0
  #DEBUG=True
fi
rm -rf $NAME/$NPROC/$SEED
mkdir -p $NAME/$NPROC/$SEED
cd $NAME/$NPROC/$SEED

if [ -z $FOAM ] ; then
  echo "OpenFOAM directory not found; assuming script is executing inside a container started by launch.sh"
else
  module load gcc/11 openmpi/gcc/4.1.6 
  source $FOAM/etc/bashrc
fi

# Set DUMP=True to save everything needed to reconstruct the linear systems
# (matrices, RHS, solutions) of every solve into the directory passed as the
# second argument to runSimulation. Requires "#define DUMP_ABSOL" in
# src/PCGBandit/PCGBandit.C. Leave empty for normal benchmarking.
DUMP=True

# Save the solver logs, every <field>-Absol dump tree, and the parallel
# cell/face/point/boundary processor addressing (needed to map the decomposed
# cells back when assembling the global matrix) into directory $1, preserving
# relative paths so Absol.py can read it directly. The bulk mesh geometry and
# field data are discarded and the run output is reset for the next call. Works
# for serial, decomposed (processorN) and collated (processorsN) runs.
dumpSave() {
  local dest=$1 p paths
  mkdir -p "$dest"
  mv log.* "$dest"/ 2>/dev/null
  mapfile -t paths < <(find . \( -type d -name '*-Absol' -prune -print \) -o \
    \( -type f -name '*ProcAddressing' -print \))
  for p in "${paths[@]}" ; do
    mkdir -p "$dest/$(dirname "$p")"
    mv "$p" "$dest/$(dirname "$p")/"
  done
  rm -rf processor[0-9]* processors[0-9]* postProcessing 2>/dev/null
  find . -mindepth 1 -maxdepth 1 -type d -name '[0-9]*' ! -name '0' ! -name '0.orig' \
    -exec rm -rf {} + 2>/dev/null
}

DEFAULT='solver PCGBandit; preconditioner separate; smootherTune yes; nCellsInCoarsestLevelTune yes; mergeLevelsTune yes; numDroptols 8; static 8;'

################################################################################
##################################   cases   ###################################
################################################################################

if [ $NAME == "boxTurb32" ] ; then

  if [ ! $DUMP ] ; then
    SWEEP=True
  fi

  cp -r $FOAM_TUTORIALS/DNS/dnsFoam/boxTurb16 boxTurb32
  cd boxTurb32

  sed -i '16a\libs ( libICTCPreconditioner.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '18a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  if [ ! $DUMP ] ; then
    sed -i 's/writeInterval   0\.25/writeInterval   100/' system/controlDict
  fi
  if [ $DEBUG ] ; then
    sed -i 's/endTime         10/endTime         0.01/' system/controlDict
  fi

  sed -i 's/16/32/g' system/blockMeshDict
  sed -i 's/nonuniform List<vector>/uniform (0 0 0);/' 0.orig/U
  sed -i '20,4119d' 0.orig/U
  sed -i 's/deltaT          0\.025/deltaT          0.005/' system/controlDict

  sed -i '21d' system/fvSolution
  sed -i 's/preconditioner  DIC;/'"$DEFAULT"'/' system/fvSolution

  runSimulation() {
    sed -i 's/'"$DEFAULT"'/'"$1"'/' system/fvSolution
    cp -r 0.orig 0
    blockMesh > log.blockMesh
    boxTurb > log.boxTurb
    echo $2 $1
    dnsFoam > log.dnsFoam
    if [ $DUMP ] ; then 
      dumpSave $2
    else 
      mv log.dnsFoam $2
      foamCleanTutorials
    fi
    sed -i 's/'"$1"'/'"$DEFAULT"'/' system/fvSolution
  }

fi

################################################################################

if [ $NAME == "pitzDaily" ] ; then

  if [ ! DUMP ] ; then
    SWEEP=True
  fi

  cp -r $FOAM_TUTORIALS/incompressible/pimpleFoam/RAS/pitzDaily .
  cd pitzDaily

  sed -i '16a\libs ( libICTCPreconditioner.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '18a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  if [ $DUMP ] ; then
    sed -i 's/writeInterval   0\.01/writeInterval   0.005/' system/controlDict
  else
    sed -i 's/writeInterval   0\.01/writeInterval   1/' system/controlDict
  fi
  if [ $DEBUG ] ; then
    sed -i 's/endTime         0\.3/endTime         0.0003/' system/controlDict
  fi

  sed -i 's/(1 3 0\.3)/(1 3 0.15)/' system/blockMeshDict
  sed -i 's/(2 4 0\.25)/(2 4 0.125)/' system/blockMeshDict
  sed -i 's/(1 1 0\.25)/(1 1 0.125)/' system/blockMeshDict
  sed -i 's/(18 30 1)/(36 60 1)/' system/blockMeshDict
  sed -i 's/(180 27 1)/(360 54 1)/' system/blockMeshDict
  sed -i 's/(180 30 1)/(360 60 1)/' system/blockMeshDict
  sed -i 's/(25 27 1)/(50 54 1)/' system/blockMeshDict
  sed -i 's/(25 30 1)/(50 60 1)/' system/blockMeshDict

  sed -i '21d' system/fvSolution
  sed -i 's/smoother         DICGaussSeidel;/'"$DEFAULT"'/' system/fvSolution

  runSimulation() {
    sed -i 's/'"$DEFAULT"'/'"$1"'/' system/fvSolution
    blockMesh > log.blockMesh
    echo $2 $1
    pimpleFoam > log.pimpleFoam
    if [ $DUMP ] ; then 
      dumpSave $2
    else 
      mv log.pimpleFoam $2
      foamCleanTutorials
    fi
    sed -i 's/'"$1"'/'"$DEFAULT"'/' system/fvSolution
  }

fi

################################################################################

if [ $NAME == "interStefanProblem" ] ; then

  mkdir interStefanProblem
  cd interStefanProblem
  for SUBFOLDER in 0.orig constant system ; do
    mkdir $SUBFOLDER
    cp -r $FOAM_TUTORIALS/verificationAndValidation/multiphase/StefanProblem/setups.orig/common/$SUBFOLDER/* $SUBFOLDER/.
    cp -r $FOAM_TUTORIALS/verificationAndValidation/multiphase/StefanProblem/setups.orig/interCondensatingEvaporatingFoam/$SUBFOLDER/* $SUBFOLDER/.
  done

  sed -i '17a\libs ( libICTCPreconditioner.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '19a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  if [ $DUMP ] ; then
    sed -i 's/writeInterval   5/writeInterval   1/' system/controlDict
  else
    sed -i 's/writeInterval   5/writeInterval   100/' system/controlDict
  fi
  if [ $DEBUG ] ; then
    sed -i 's/endTime         50/endTime         1.4/' system/controlDict
  fi

  sed -i 's/(400 2 1)/(800 4 1)/' system/blockMeshDict

  sed -i '45d' system/fvSolution
  sed -i '/solver *PCG;/d' system/fvSolution
  sed -i '/smoother *DIC;/d' system/fvSolution
  sed -i 's/preconditioner *DIC;/'"$DEFAULT"'/' system/fvSolution

  runSimulation() {
    sed -i 's/'"$DEFAULT"'/'"$1"'/' system/fvSolution
    cp -r 0.orig 1.36
    blockMesh > log.blockMesh
    renumberMesh -overwrite -constant > log.renumberMesh
    checkMesh > log.checkMesh
    setAlphaField > log.setAlphaField
    echo $2 $1
    if (( $NPROC > 1 )) ; then
      sed -i 's/numberOfSubdomains.*;/numberOfSubdomains '"$NPROC"';/' system/decomposeParDict
      decomposePar > log.decomposePar
      if [ $DEBUG ] || [ -z $FOAM ] ; then
        mpirun -n $NPROC interCondensatingEvaporatingFoam -parallel > log.interCondensatingEvaporatingFoam
      else
        srun -n $NPROC interCondensatingEvaporatingFoam -parallel > log.interCondensatingEvaporatingFoam
      fi
    else
      interCondensatingEvaporatingFoam > log.interCondensatingEvaporatingFoam
    fi
    if [ $DUMP ] ; then 
      dumpSave $2
    else 
      mv log.interCondensatingEvaporatingFoam $2
      foamCleanTutorials
    fi
    sed -i 's/'"$1"'/'"$DEFAULT"'/' system/fvSolution
  }

fi

################################################################################

if [ $NAME == "porousDamBreak" ] ; then

  if [ $NPROC -ge 16 ] && [ ! $DUMP ] ; then
    SWEEP=True
  fi

  cp -r $FOAM_TUTORIALS/verificationAndValidation/multiphase/interIsoFoam/porousDamBreak porousDamBreak
  cd porousDamBreak

  sed -i '16a\libs ( libICTCPreconditioner.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '18a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  sed -i 's/writeInterval   0\.05/writeInterval   0.01/' system/controlDict
  if [ ! $DUMP ] ; then
    sed -i 's/writeInterval   0\.05/writeInterval   10/' system/controlDict
  fi
  if [ $DEBUG ] ; then
    sed -i 's/endTime         4/endTime         0.004/' system/controlDict
  fi

  sed -i 's/(75 93 1)/(150 186 1)/' system/blockMeshDict
  sed -i 's/(73 93 1)/(146 186 1)/' system/blockMeshDict
  sed -i 's/(76 93 1)/(152 186 1)/' system/blockMeshDict
  sed -i 's/(75 53 1)/(150 106 1)/' system/blockMeshDict
  sed -i 's/(73 53 1)/(146 106 1)/' system/blockMeshDict
  sed -i 's/(76 53 1)/(152 106 1)/' system/blockMeshDict

  sed -i '/solver *PCG;/d' system/fvSolution
  sed -i 's/preconditioner *DIC;/'"$DEFAULT"'/' system/fvSolution

  runSimulation() {
    sed -i 's/'"$DEFAULT"'/'"$1"'/' system/fvSolution
    cp -r 0.orig 0
    blockMesh > log.blockMesh
    setFields > log.setFields
    echo $2 $1
    if (( $NPROC > 1 )) ; then
      sed -i 's/numberOfSubdomains.*;/numberOfSubdomains '"$NPROC"';/' system/decomposeParDict
      for n in $(seq 0 $NPROC) ; do
        if (( $((n * n)) == $NPROC )) ; then
          sed -i 's/n *(.*1);/n ('"$n"' '"$n"' 1);/' system/decomposeParDict
          break
        fi
      done
      decomposePar > log.decomposePar
      if [ $DEBUG ] || [ -z $FOAM ] ; then
        mpirun -n $NPROC interIsoFoam -parallel > log.interIsoFoam
      else
        srun -n $NPROC interIsoFoam -parallel > log.interIsoFoam
      fi
    else
      interIsoFoam > log.interIsoFoam
    fi
    if [ $DUMP ] ; then 
      dumpSave $2
    else 
      mv log.interIsoFoam $2
      foamCleanTutorials
    fi
    sed -i 's/'"$1"'/'"$DEFAULT"'/' system/fvSolution
  }

fi

################################################################################

if [ $NAME == "closedPipe" ] ; then

  echo "running FreeMHD simulation; assuming examples/FreeMHD.zip has been unzipped"
  wmake ../../../FreeMHD/solvers/epotMultiRegionInterFoam

  cp -r ../../../FreeMHD/closedPipe .
  cd closedPipe
  wmake libso dynamicCode/outletUxB

  sed -i '19a\libs ( libICTCPreconditioner.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '21a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  if [ $DUMP ] ; then
    sed -i 's/writeInterval   1e-5/writeInterval   5e-4/' system/controlDict
  else
    sed -i 's/writeInterval   1e-5/writeInterval   1e-1/' system/controlDict
  fi
  if [ $DEBUG ] ; then
    sed -i 's/endTime         0.025/endTime         2.5e-6/' system/controlDict
  fi

  for region in $(foamListRegions) ; do
    sed -i '/solver *PCG;/d' system/$region/fvSolution
    sed -i 's/preconditioner *DIC;/'"$DEFAULT"'/' system/$region/fvSolution
  done

  runSimulation() {
    cp -r 0.orig 0
    blockMesh > log.blockMesh
    topoSet > log.topoSet
    splitMeshRegions -cellZonesOnly -overwrite -fileHandler collated > log.splitMeshRegions
    for region in $(foamListRegions) ; do
      sed -i 's/'"$DEFAULT"'/'"$1"'/' system/$region/fvSolution
      changeDictionary -region $region -fileHandler collated > log.changeDictionary.$region
      setExprFields -region $region -fileHandler collated > log.setExprFields.$region
    done

    echo $2 $1
    if (( $NPROC > 1 )) ; then
      sed -i 's/numberOfSubdomains.*;/numberOfSubdomains '"$NPROC"';/' system/decomposeParDict
      for region in $(foamListRegions) ; do
        sed -i 's/numberOfSubdomains.*;/numberOfSubdomains '"$NPROC"';/' system/$region/decomposeParDict
      done
      decomposePar -allRegions -force -fileHandler collated > log.decomposePar
      if [ $DEBUG ] || [ -z $FOAM ] ; then
        mpirun -n $NPROC epotMultiRegionInterFoam -parallel > log.epotMultiRegionInterFoam
      else
        srun -n $NPROC epotMultiRegionInterFoam -parallel > log.epotMultiRegionInterFoam
      fi
    else
      epotMultiRegionInterFoam > log.epotMultiRegionInterFoam
    fi
    if [ $DUMP ] ; then
      dumpSave $2
    else
      mv log.epotMultiRegionInterFoam $2
      rm -rf log.* 0 constant/cellToRegion postProcessing
      for region in $(foamListRegions) ; do rm -rf constant/$region/polyMesh ; done
    fi
    for region in $(foamListRegions) ; do
      sed -i 's/'"$1"'/'"$DEFAULT"'/' system/$region/fvSolution
    done
  }

fi

################################################################################

if [ $NAME == 'fringingBField' ] ; then

  echo "running FreeMHD simulation; assuming examples/FreeMHD.zip has been unzipped"
  wmake ../../../FreeMHD/solvers/epotMultiRegionInterFoam

  cp -r ../../../FreeMHD/fringingBField .
  cd fringingBField
  wmake libso dynamicCode/outletUxB
  
  sed -i '19a\libs ( libICTCPreconditioner.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '21a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  if [ ! $DUMP ] ; then
    sed -i 's/writeInterval   0.1/writeInterval   1/' system/controlDict
  fi
  if [ $DEBUG ] ; then
    sed -i 's/endTime         0.5/endTime         5e-6/' system/controlDict
  fi

  for region in $(foamListRegions) ; do
    sed -i '/solver *PCG;/d' system/$region/fvSolution
    sed -i 's/preconditioner *DIC;/'"$DEFAULT"'/' system/$region/fvSolution
  done

  runSimulation() {
    cp -r 0.orig 0
    blockMesh > log.blockMesh
    splitMeshRegions -cellZonesOnly -overwrite -fileHandler collated > log.splitMeshRegions
    for region in $(foamListRegions) ; do
      sed -i 's/'"$DEFAULT"'/'"$1"'/' system/$region/fvSolution
      changeDictionary -region $region -fileHandler collated > log.changeDictionary.$region
      setExprFields -region $region -fileHandler collated > log.setExprFields.$region
    done

    echo $2 $1
    if (( $NPROC > 1 )) ; then
      sed -i 's/numberOfSubdomains.*;/numberOfSubdomains '"$NPROC"';/' system/decomposeParDict
      for region in $(foamListRegions) ; do
        sed -i 's/numberOfSubdomains.*;/numberOfSubdomains '"$NPROC"';/' system/$region/decomposeParDict
      done
      decomposePar -allRegions -force -fileHandler collated > log.decomposePar
      if [ $DEBUG ] || [ -z $FOAM ] ; then
        mpirun -n $NPROC epotMultiRegionInterFoam -parallel > log.epotMultiRegionInterFoam
      else
        srun -n $NPROC epotMultiRegionInterFoam -parallel > log.epotMultiRegionInterFoam
      fi
    else
      epotMultiRegionInterFoam > log.epotMultiRegionInterFoam
    fi
    if [ $DUMP ] ; then
      dumpSave $2
    else
      mv log.epotMultiRegionInterFoam $2
      rm -rf log.* 0 constant/cellToRegion postProcessing
      for region in $(foamListRegions) ; do rm -rf constant/$region/polyMesh ; done
    fi
    for region in $(foamListRegions) ; do
      sed -i 's/'"$1"'/'"$DEFAULT"'/' system/$region/fvSolution
    done
  }

fi

# controlDict is stable after the case setup above, so apply the dump settings
# once here (cwd is the case directory) to cover every runSimulation call.
if [ $DUMP ] ; then
  sed -i 's/^\( *purgeWrite *\).*/\10;/' system/controlDict
  sed -i 's/^\( *writeFormat *\).*/\1binary;/' system/controlDict
fi

################################################################################
###################################   runs   ###################################
################################################################################

if [ $SWEEP ] ; then
  STATIC=`seq 0 32`   # try all configurations
else
  STATIC='8 21'       # try DIC & GAMG(DIC+GS)
fi

for S in $STATIC ; do
  LOGFILE=../staticPCG_"$S"
  SOLVER='solver PCGBandit; preconditioner separate; smootherTune yes; nCellsInCoarsestLevelTune yes; mergeLevelsTune yes; numDroptols 8; static '"$S"'; backstop 10000; cacheAgglomeration no;'
  runSimulation "$SOLVER" $LOGFILE
done

################################################################################

SOLVER='solver PCGBandit; preconditioner separate; smootherTune yes; nCellsInCoarsestLevelTune yes; mergeLevelsTune yes; numDroptols 8; backstop 10000; cacheAgglomeration no;'
runSimulation "$SOLVER" ../PCGBandit

################################################################################

cd ..
rm -rf $NAME
