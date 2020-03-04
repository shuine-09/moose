#!/bin/bash

#
# Take multiple *.png sequences outputted by paraview and generate a
# single side by side gif animation
# Arguments: base names of the image series
#

# output file name
outname=`echo $@ | sed 's/ /_/g'`.gif

# montage the series together
mkdir .tmp
for file in $1.*.png
do
  # get number
  num=`echo $file | cut -d. -f2`
  arg=''
  for base in $@
  do
    arg=$arg' '${base}.${num}.png
  done
  convert +append $arg .tmp/${num}.gif
done

# animate
for anim in `ls | cut -d. -f1|sort -u`
do
  gifsicle .tmp/*.gif --colors 256 -l > $outname
done

# cleanup
rm -rf .tmp
