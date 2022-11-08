#!/bin/bash -e

chop() {
  swc modify $1 -i $(swc find $1 -p 3 -g 0) -u -o $2
}

name=$(basename $1 .swc)
cp $1 $name-00.swc
echo -n chopping ... ''
for i in $(seq 0 29); do
  echo -n $((i+1)) ''
  i0=$(printf '%02d' $i)
  i1=$(printf '%02d' $((i+1)))
  chop $name-$i0.swc $name-$i1.swc
done
echo done
