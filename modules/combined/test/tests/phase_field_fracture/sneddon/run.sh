#!/bin/bash
#a=(0.1 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.45 0.5 0.55 0.6)
a=(0.4 0.45 0.5 0.55 0.6)
outdirname='critical_p'
#\rm -r $outdirname
#mkdir $outdirname

for i in ${!a[@]};
do
  mpiexec -n 8 ../../../../combined-opt -i sneddon_disp.i a=${a[$i]} damage:a=${a[$i]} Outputs/file_base=$outdirname/${a[$i]}_out 2>&1 | tee $outdirname/${a[$i]}.log
done
