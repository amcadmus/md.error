#!/bin/bash

target=$*

for i in $target;
do
    echo "# processing $i"
    p1=`echo $i | cut -d '.' -f 1`
    p2=`echo $i | cut -d '.' -f 2`
    p3=`echo $i | cut -d '.' -f 3`
    out=$p1.$p2.$p3.gif
    sed -e "/^set out/s/'.*'/'$out'/g" error.one.gp |\
    sed -e "/^spl/s/'.*'/'$i'/g" > tmp.gp
    gnuplot tmp.gp
done

gifsicle --delay=10 $p1.*.*.gif > error.out.gif
