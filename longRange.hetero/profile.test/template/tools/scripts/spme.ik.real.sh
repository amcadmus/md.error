#!/bin/bash

source parameters.sh

if test ! -d $record_dir; then
    echo "no record dir $record_dir"
    exit
fi
if test ! -d $errors_dir; then
    echo "# no errors dir $errors_dir, make it"
    mkdir -p $errors_dir
fi

mylog=spme.ik.real.log
rm -f $mylog
spme_ik_record=spme.ik.record
rm -f $spme_ik_record

make -C ./tools/real.error/ &>> make.log

for i in `ls $record_dir | grep conf | grep gro`;
do
    iname=`echo $i | cut -d '.' -f 2`
    echo "# #############################################################################" &>> $mylog
    echo "# doing $iname" &>> $mylog
    echo "# #############################################################################" &>> $mylog
    ./tools/real.error/spme --config $record_dir/$i --method bspline --beta $beta --rcut $cal_rcut --grid-number $cal_K --order $cal_order --output-force force.spme.ik.out &>> $mylog
    
    tools/real.error/compare.force --config $record_dir/$i --refh $real_h --force-0 $record_dir/force.ref.$iname.out --force-1 force.spme.ik.out -o error.spme.ik.$iname.out -m meanf.spme.ik.$iname.out &>> $mylog

    mv error.spme.ik.$iname.out meanf.spme.ik.$iname.out $errors_dir
    echo "$errors_dir/error.spme.ik.$iname.out $errors_dir/meanf.spme.ik.$iname.out" >> $spme_ik_record
done
mv -f $spme_ik_record $errors_dir

echo "# #############################################################################" &>> $mylog
echo "# combine results" &>> $mylog
echo "# #############################################################################" &>> $mylog
./tools/real.error/combine.results --record-file $errors_dir/$spme_ik_record -o $errors_dir/real.spme.ik.error.out -m $errors_dir/real.spme.ik.meanf.out &>> $mylog

rm -f force.spme.ik.out
mv -f make.log $mylog $errors_dir
