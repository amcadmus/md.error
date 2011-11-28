#!/bin/bash

source parameters.sh

mylog=ewald.log
rm -f $mylog
ewald_record=ewald.record
rm -f $ewald_record

make -C ./tools/real.error/ &>> make.log

for i in `ls $record_dir | grep conf | grep gro`;
do
    iname=`echo $i | cut -d '.' -f 2`
    echo "# #############################################################################" &>> $mylog
    echo "# doing $iname" &>> $mylog
    echo "# #############################################################################" &>> $mylog
    ./tools/real.error/ewaldSum --config $record_dir/$i --beta $beta --rcut $cal_rcut --grid-number $cal_K --output-force force.ewald.out &>> $mylog
    
    tools/real.error/compare.force --config $record_dir/$i --force-0 $record_dir/force.ref.$iname.out --force-1 force.ewald.out -o error.ewald.$iname.out -m meanf.ewald.$iname.out

    echo "error.ewald.$iname.out meanf.ewald.$iname.out" >> $ewald_record
    mv error.ewald.$iname.out meanf.ewald.$iname.out $record_dir
done

rm -f force.ewald.out
mv $ewald_record $record_dir

