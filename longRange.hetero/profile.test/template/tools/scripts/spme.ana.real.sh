#!/bin/bash

source parameters.sh

mylog=spme.ana.real.log
rm -f $mylog
spme_ana_record=spme.ana.record
rm -f $spme_ana_record

make -C ./tools/real.error/ &>> make.log

for i in `ls $record_dir | grep conf | grep gro`;
do
    iname=`echo $i | cut -d '.' -f 2`
    echo "# #############################################################################" &>> $mylog
    echo "# doing $iname" &>> $mylog
    echo "# #############################################################################" &>> $mylog
    ./tools/real.error/spme --config $record_dir/$i --method fbspline --beta $beta --rcut $cal_rcut --grid-number $cal_K --order $cal_order --output-force force.spme.ana.out &>> $mylog
    
    tools/real.error/compare.force --config $record_dir/$i --refh $real_h --force-0 $record_dir/force.ref.$iname.out --force-1 force.spme.ana.out -o error.spme.ana.$iname.out -m meanf.spme.ana.$iname.out &>> $mylog

    mv error.spme.ana.$iname.out meanf.spme.ana.$iname.out $record_dir
    echo "$record_dir/error.spme.ana.$iname.out $record_dir/meanf.spme.ana.$iname.out" >> $spme_ana_record
done
mv -f $spme_ana_record $record_dir

echo "# #############################################################################" &>> $mylog
echo "# combine results" &>> $mylog
echo "# #############################################################################" &>> $mylog
./tools/real.error/combine.results --record-file $record_dir/$spme_ana_record -o $record_dir/real.spme.ana.error.out -m $record_dir/real.spme.ana.meanf.out &>> $mylog

# # esti the dir error
# ./tools/analyze/error.dir -t $record_dir/traj.xtc --refh $real_h --beta $beta --rcut $cal_rcut
# mv -f rho.x.avg.out $record_dir/
# mv -f error.out $record_dir/esti.dir.error.out
# mv -f meanf.out $record_dir/esti.dir.meanf.out

rm -f force.spme.ana.out
mv -f $mylog $record_dir
