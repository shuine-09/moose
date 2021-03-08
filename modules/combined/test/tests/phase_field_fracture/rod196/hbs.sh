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
mpiexec $HOME/projects/moose5_lemhi/modules/combined/combined-opt -i bubble_disp.i AuxVariables/c/InitialCondition/radius=0.25  damage:AuxVariables/c/InitialCondition/radius=0.25  Functions/pressure/data_file='bubble_pressure_r0.25_gas100_por10_ext0_rod196.csv' damage:Functions/pressure/data_file='bubble_pressure_r0.25_gas100_por10_ext0_rod196.csv' BCs/Pressure/coolantPressure/factor=0  Outputs/file_base='bubble_pressure_r0.25_gas100_por10_ext0_rod196_out' >> bubble_pressure_r0.25_gas100_por10_ext0_rod196
date >> hbs
