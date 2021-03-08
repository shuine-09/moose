#!/bin/bash
#a=(20 10 5 2.5 1.25 0.625)
#b=(200 400 800 1600 3200 6400)
#c=(1 2 4 8 16 32 64)
a=(1.25 0.625)
b=(3200 6400)
c=(32 64)
outdirname='output2'
#\rm -r $outdirname
#mkdir $outdirname

for i in ${!a[@]};
do
  mpiexec -n 20 ../../../../combined-opt -i bar_disp.i Functions/fracture_pressure/value=0.4 damage:Functions/fracture_pressure/value=0.4 Mesh/gen/nx=${b[$i]} Mesh/gen/ny=${c[$i]} damage:Mesh/gen/nx=${b[$i]} damage:Mesh/gen/ny=${c[$i]} Materials/pfbulkmat2/prop_values=${a[$i]} damage:Materials/pfbulkmat2/prop_values=${a[$i]}  Outputs/file_base=$outdirname/${b[$i]}_out 2>&1 | tee $outdirname/${b[$i]}.log
done
