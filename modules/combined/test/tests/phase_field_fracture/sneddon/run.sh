#!/bin/bash
a=(0.20 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4)
outdirname='critical_p'
\rm -r $outdirname
mkdir $outdirname

for i in ${!a[@]};
do
  mpiexec -n 8 ../../../../combined-opt -i sneddon_disp.i a=${a[$i]} damage:a=${a[$i]} Outputs/file_base=$outdirname/${a[$i]}_out 2>&1 | tee $outdirname/${a[$i]}.log
done
