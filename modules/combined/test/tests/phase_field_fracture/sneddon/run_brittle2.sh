#!/bin/bash
a=(0.1 0.15 0.20 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6)
outdirname='critical_l004'
\rm -r $outdirname
mkdir $outdirname

for i in ${!a[@]};
do
  mpiexec -n 10 ../../../../combined-opt -i sneddon_disp.i a=${a[$i]} damage:a=${a[$i]} Materials/pfbulkmat2/prop_values=0.04 damage:Materials/pfbulkmat2/prop_values=0.04 Outputs/file_base=$outdirname/${a[$i]}_out 2>&1 | tee $outdirname/${a[$i]}.log
done
