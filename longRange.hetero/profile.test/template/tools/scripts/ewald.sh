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
    
    tools/real.error/compare.force --config $record_dir/$i --refh $real_h --force-0 $record_dir/force.ref.$iname.out --force-1 force.ewald.out -o error.ewald.$iname.out -m meanf.ewald.$iname.out &>> $mylog

    mv error.ewald.$iname.out meanf.ewald.$iname.out $record_dir
    echo "$record_dir/error.ewald.$iname.out $record_dir/meanf.ewald.$iname.out" >> $ewald_record
done
mv -f $ewald_record $record_dir

echo "# #############################################################################" &>> $mylog
echo "# combine results" &>> $mylog
echo "# #############################################################################" &>> $mylog
./tools/real.error/combine.results --record-file $record_dir/$ewald_record -o $record_dir/real.error.out -m $record_dir/real.meanf.out &>> $mylog

# # esti the dir error
# ./tools/analyze/error.dir -t $record_dir/traj.xtc --refh $real_h --beta $beta --rcut $cal_rcut
# mv -f rho.x.avg.out $record_dir/
# mv -f error.out $record_dir/esti.dir.error.out
# mv -f meanf.out $record_dir/esti.dir.meanf.out

rm -f force.ewald.out
mv -f $mylog $record_dir
