#!/bin/bash
a=(0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3)
outdirname='critical_l002'
\rm -r $outdirname
mkdir $outdirname

for i in ${!a[@]};
do
  mpiexec -n 16 ../../../../combined-opt -i sneddon_disp.i a=${a[$i]} damage:a=${a[$i]} Materials/pfbulkmat2/prop_values=0.02 damage:Materials/pfbulkmat2/prop_values=0.02 Outputs/file_base=$outdirname/${a[$i]}_out 2>&1 | tee $outdirname/${a[$i]}.log
done
