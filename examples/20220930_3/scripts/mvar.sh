#!/bin/bash -e
for f in *-var0.swc; do
  n=`basename $f -var0.swc`
  echo -n $n ''
  for i in `seq 1 2 7`; do
      s=`echo 1.05 + 0.05*$i | bc`
      echo -n $s ''
      swc modify $f -w360 -p3
      swc modify mod.swc -s $s $s $s
      swc repair mod.swc -r 3 -n -o $n-var$i.swc
      rm mod.swc
  done
  for i in `seq 2 2 8`; do
      s=`echo 1.0 - 0.05*$i | bc`
      echo -n $s ''
      swc modify $f -w360 -p3
      swc modify mod.swc -s $s $s $s
      swc repair mod.swc -r 3 -n -o $n-var$i.swc
      rm mod.swc
  done
  echo
done
