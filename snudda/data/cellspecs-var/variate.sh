#!/bin/bash
for d in dspn fs ispn; do
  echo $d
  pushd $d > /dev/null
    for c in *; do 
      f=`basename $c/*.swc`
      n=`basename $f .swc`
      echo $c
      for i in {0..8}; do 
        echo -n $i ''
        cp -r $c $c-var$i
        s=`echo 0.8+0.05*$i | bc -l`
        mmod -c $s,$s,$s -o tmp.swc $c/$f
        mmod -p3 -i20 -o $c-var$i/$n-var$i.swc tmp.swc
        rm $c-var$i/$f
      done
      rm -r $c tmp.swc 
      echo
    done
  popd > /dev/null
  echo
done
