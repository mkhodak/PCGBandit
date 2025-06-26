#!/bin/bash
#SBATCH --job-name=pitzDaily
#SBATCH --array=[0-2]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --nodelist=della-r3c[1-4]n[1-16]
#SBATCH --nodes=1

hostname
if [ -z $FOAM ] ; then
  echo "OpenFOAM directory not found"
fi

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

module load gcc/11 openmpi/gcc/4.1.6 
source $FOAM/etc/bashrc

BACKSTOP=' maxIter 2000; cacheAgglomeration no;'
DEFAULT='solver PCGBandit; preconditioner separate; smootherTune yes; nCellsInCoarsestLevelTune yes; mergeLevelsTune yes; numDroptols 8; static 8;'

################################################################################
##################################   cases   ###################################
################################################################################

if [ $NAME == "pitzDaily" ] ; then

  SWEEP=True

  cp -r $FOAM_TUTORIALS/incompressible/pimpleFoam/RAS/pitzDaily .
  cd pitzDaily

  sed -i '16a\libs ( libICTCPreconditioner.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '18a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  sed -i 's/writeInterval   0\.01/writeInterval   1/' system/controlDict
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
    mv log.pimpleFoam $2
    foamCleanTutorials
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
  sed -i 's/writeInterval   5/writeInterval   100/' system/controlDict
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
      if [ $DEBUG ] ; then
        mpirun -n $NPROC interCondensatingEvaporatingFoam -parallel > log.interCondensatingEvaporatingFoam
      else
        srun -n $NPROC interCondensatingEvaporatingFoam -parallel > log.interCondensatingEvaporatingFoam
      fi
    else
      interCondensatingEvaporatingFoam > log.interCondensatingEvaporatingFoam
    fi
    mv log.interCondensatingEvaporatingFoam $2
    foamCleanTutorials
    sed -i 's/'"$1"'/'"$DEFAULT"'/' system/fvSolution
  }

fi

################################################################################

if [ $NAME == "closedPipe" ] ; then

  cp -r ../../../FreeMHD/closedPipe .
  cd closedPipe
  wmake libso dynamicCode/outletUxB

  sed -i '19a\libs ( libICTCPreconditioner.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '21a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  sed -i 's/writeInterval   1e-5/writeInterval   1e-1/' system/controlDict
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
      if [ $DEBUG ] ; then
        mpirun -n $NPROC epotMultiRegionInterFoam -parallel > log.epotMultiRegionInterFoam
      else
        srun -n $NPROC epotMultiRegionInterFoam -parallel > log.epotMultiRegionInterFoam
      fi
    else
      epotMultiRegionInterFoam > log.epotMultiRegionInterFoam
    fi
    mv log.epotMultiRegionInterFoam $2

    rm -rf log.*
    rm -rf 0
    rm -rf constant/cellToRegion
    rm -rf postProcessing
    for region in $(foamListRegions) ; do
      rm -rf constant/$region/polyMesh
      sed -i 's/'"$1"'/'"$DEFAULT"'/' system/$region/fvSolution
    done
  }

fi

################################################################################

if [ $NAME == 'fringingBField' ] ; then

  cp -r ../../../FreeMHD/fringingBField .
  cd fringingBField
  wmake libso dynamicCode/outletUxB
  
  sed -i '19a\libs ( libICTCPreconditioner.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '21a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  sed -i 's/writeInterval   0.1/writeInterval   1/' system/controlDict
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
      if [ $DEBUG ] ; then
        mpirun -n $NPROC epotMultiRegionInterFoam -parallel > log.epotMultiRegionInterFoam
      else
        srun -n $NPROC epotMultiRegionInterFoam -parallel > log.epotMultiRegionInterFoam
      fi
    else
      epotMultiRegionInterFoam > log.epotMultiRegionInterFoam
    fi
    mv log.epotMultiRegionInterFoam $2

    rm -rf log.*
    rm -rf 0
    rm -rf constant/cellToRegion
    rm -rf postProcessing
    for region in $(foamListRegions) ; do
      rm -rf constant/$region/polyMesh
      sed -i 's/'"$1"'/'"$DEFAULT"'/' system/$region/fvSolution
    done
  }

fi

################################################################################
###################################   runs   ###################################
################################################################################

if [ $SWEEP ] ; then
  STATIC=`seq 0 32`
else
  STATIC='8 21'
fi

for S in $STATIC ; do
  LOGFILE=../staticPCG_"$S"
  SOLVER='solver PCGBandit; preconditioner separate; smootherTune yes; nCellsInCoarsestLevelTune yes; mergeLevelsTune yes; numDroptols 8; static '"$S"';'"$BACKSTOP"
  runSimulation "$SOLVER" $LOGFILE
  echo "$SOLVER" >> $LOGFILE
done

################################################################################

SOLVER='solver PCGBandit; preconditioner separate; smootherTune yes; nCellsInCoarsestLevelTune yes; mergeLevelsTune yes; numDroptols 8;'"$BACKSTOP"
runSimulation "$SOLVER" ../PCGBandit

################################################################################

cd ..
rm -rf $NAME
