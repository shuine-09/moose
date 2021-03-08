#!/bin/bash
a=(0.0 0.2 0.4 0.6 0.8 1.0)
outdirname='output_new'
\rm -r $outdirname
mkdir $outdirname

for i in ${!a[@]};
do
  mpiexec -n 4 ../../../../combined-opt -i bar_disp.i Functions/fracture_pressure/value=${a[$i]} damage:Functions/fracture_pressure/value=${a[$i]} Outputs/file_base=$outdirname/${a[$i]}_out 2>&1 | tee $outdirname/${a[$i]}.log
done
