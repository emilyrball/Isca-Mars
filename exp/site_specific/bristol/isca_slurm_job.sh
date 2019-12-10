#!/bin/sh

#PBS -l select=1:ncpus=16
#PBS -l walltime=1:00:00
#PBS -o /home/xz19136/Isca_jobs/test.o
#PBS -e /home/xz19136/Isca_jobs/test.e

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

module purge
source $HOME/.bashrc
module load lang/python/anaconda/3.7-2019.03.biopython
source $GFDL_BASE/src/extra/env/bluepebble
source activate isca_env

$HOME/.conda/envs/isca_env/bin/python $GFDL_BASE/exp/test_cases/held_suarez/held_suarez_test_case.py
