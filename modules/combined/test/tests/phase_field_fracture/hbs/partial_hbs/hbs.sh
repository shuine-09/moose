#!/bin/bash
#PBS -m abe
#PBS -N hbs
#PBS -l select=4:ncpus=40:mpiprocs=20
#PBS -l walltime=72:00:00
#PBS -P neams

JOB_NUM=${PBS_JOBID%%\.*}

cd $PBS_O_WORKDIR

module purge
module load pbs
module load use.moose PETSc

\rm -f hbs
date > hbs
mpiexec $HOME/projects/moose5_lemhi/modules/combined/combined-opt -i test_disp.i  UserObjects/soln/mesh=2020_12_07_HBS_bub_vol_frac1.e  damage:UserObjects/soln/mesh=2020_12_07_HBS_bub_vol_frac1.e Outputs/file_base='test1'  >> test1
date >> hbs
