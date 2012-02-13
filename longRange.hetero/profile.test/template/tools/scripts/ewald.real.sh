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
rm -f $errors_dir/parameters.ewald.real.sh
cp parameters.sh $errors_dir/parameters.ewald.real.sh

mylog=ewald.real.log
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
    if echo $project_name | grep one_peak &> /dev/null; then
	./tools/real.error/ewaldSum --config $record_dir/$i --charge-table $record_dir/charge.tab --beta $beta --rcut $cal_rcut --kx $cal_Kx --ky $cal_Ky --kz $cal_Kz --output-force force.ewald.out &>> $mylog
    else
	./tools/real.error/ewaldSum --config $record_dir/$i --beta $beta --rcut $cal_rcut --kx $cal_Kx --ky $cal_Ky --kz $cal_Kz --output-force force.ewald.out &>> $mylog
    fi
    
    tools/real.error/compare.force --config $record_dir/$i --refh $real_h --force-0 $record_dir/force.ref.$iname.out --force-1 force.ewald.out -o error.ewald.$iname.out -m meanf.ewald.$iname.out &>> $mylog

    mv error.ewald.$iname.out meanf.ewald.$iname.out $errors_dir
    echo "$errors_dir/error.ewald.$iname.out $errors_dir/meanf.ewald.$iname.out" >> $ewald_record
done
mv -f $ewald_record $errors_dir

echo "# #############################################################################" &>> $mylog
echo "# combine results" &>> $mylog
echo "# #############################################################################" &>> $mylog
./tools/real.error/combine.results --record-file $errors_dir/$ewald_record -o $errors_dir/real.ewald.error.out -m $errors_dir/real.ewald.meanf.out &>> $mylog

# # esti the dir error
# ./tools/analyze/error.dir -t $record_dir/traj.xtc --refh $real_h --beta $beta --rcut $cal_rcut
# mv -f rho.x.avg.out $record_dir/
# mv -f error.out $record_dir/esti.dir.error.out
# mv -f meanf.out $record_dir/esti.dir.meanf.out

rm -f force.ewald.out
mv -f make.log $mylog $errors_dir
