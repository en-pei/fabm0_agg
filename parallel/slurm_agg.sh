#!/bin/tcsh
#SBATCH --job-name=testagg
#SBATCH --ntasks=10
#SBATCH --time=05:00:00
#SBATCH --mail-user=enpei.li@hereon.de
#SBATCH --mail-type=ALL
#SBATCH --account=KSE
#SBATCH --output=job.o%j
#SBATCH --error=job.e%j
#SBATCH --array=0-2
#SBATCH --partition=pAll,pCluster

module load applications/utils/cmake3-20.0
module load compilers/intel/2020.1.217
module load intelmpi
module load hdf5/1.10.5
module load netcdf/4.7.0
module load pnetcdf
module load applications/python/3.8

/gpfs/work/lie/bale2002/test.sh 10 $SLURM_ARRAY_TASK_ID
