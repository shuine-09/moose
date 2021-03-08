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
mpiexec $HOME/projects/moose5_lemhi/modules/combined/combined-opt -i bubble_disp.i Functions/pressure/value=0 damage:Functions/pressure/value=0 Outputs/file_base='pressure_0_out' >> pressure_0
date >> hbs
