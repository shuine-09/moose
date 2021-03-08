#!/bin/bash
#a=(0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3)
a=(0.225 0.3)
outdirname='critical_p_0p2_cz_2'
#\rm -r $outdirname
#mkdir $outdirname

for i in ${!a[@]};
do
  mpiexec -n 8 ../../../../combined-opt -i sneddon_disp_cz.i a=${a[$i]} damage:a=${a[$i]} Materials/degradation/constant_expressions='1 0.2 1e-8' damage:Materials/degradation/constant_expressions='1 0.2 1e-8' Outputs/file_base=$outdirname/${a[$i]}_out 2>&1 | tee $outdirname/${a[$i]}.log
done
