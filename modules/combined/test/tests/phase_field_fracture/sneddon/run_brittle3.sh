#!/bin/bash
#a=(0.05 0.0625 0.075 0.0875 0.1 0.1125 0.125 0.1375 0.15)
a=(0.0375)
outdirname='critical_l001'
#\rm -r $outdirname
#mkdir $outdirname

for i in ${!a[@]};
do
  mpiexec -n 8 ../../../../combined-opt -i sneddon_disp.i a=${a[$i]} damage:a=${a[$i]} Adaptivity/max_h_level=6 Adaptivity/initial_steps=6 Materials/pfbulkmat2/prop_values=0.01 damage:Materials/pfbulkmat2/prop_values=0.01 damage:Adaptivity/max_h_level=6 damage:Adaptivity/initial_steps=6 Outputs/file_base=$outdirname/${a[$i]}_out 2>&1 | tee $outdirname/${a[$i]}.log
done
