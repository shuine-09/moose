#!/bin/bash
a=(0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3)
outdirname='critical_p_0p15_cz_l001'
#\rm -r $outdirname
#mkdir $outdirname

for i in ${!a[@]};
do
  mpiexec -n 12 ../../../../combined-opt -i sneddon_disp_cz.i a=${a[$i]} damage:a=${a[$i]} Materials/degradation/constant_expressions='1 0.15 1e-8' damage:Materials/degradation/constant_expressions='1 0.15 1e-8' Adaptivity/max_h_level=6 Adaptivity/initial_steps=6 Materials/pfbulkmat/prop_values='0.001 0.01 1e-6' damage:Materials/pfbulkmat/prop_values='0.001 0.01 1e-6' damage:Adaptivity/max_h_level=6 damage:Adaptivity/initial_steps=6 Outputs/file_base=$outdirname/${a[$i]}_out 2>&1 | tee $outdirname/${a[$i]}.log
done
