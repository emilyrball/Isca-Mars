#!/bin/bash -l

#SBATCH --job-name=compile
#SBATCH --partition=cpu
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#number of tasks ~ processes per node
#SBATCH --ntasks-per-node=1
#number of cpus (cores) per task (process)
#SBATCH --cpus-per-task=1
#SBATCH --output=compile_%j.o

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

module purge
source $HOME/.bashrc
module load ifort
source $GFDL_BASE/src/extra/env/bristol-bc4
source activate isca_env

cd $GFDL_BASE/postprocessing/plevel_interpolation
./compile_plev_interpolation.sh
# $HOME/.conda/envs/isca_env/bin/python $GFDL_BASE/exp/socrates_mars/socrates_mars.py
