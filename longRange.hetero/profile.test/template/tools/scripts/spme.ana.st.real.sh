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
rm -f $errors_dir/parameters.spme.ana.st.real.sh
cp parameters.sh $errors_dir/parameters.spme.ana.st.real.sh

mylog=spme.ana.st.real.log
rm -f $mylog
spme_ana_st_record=spme.ana.st.record
rm -f $spme_ana_st_record

make -j8 -C ./tools/real.error/ &>> make.log

for i in `ls $record_dir | grep conf | grep gro`;
do
    iname=`echo $i | cut -d '.' -f 2`
    echo "# #############################################################################" &>> $mylog
    echo "# doing $iname" &>> $mylog
    echo "# #############################################################################" &>> $mylog
    ./tools/real.error/spme --config $record_dir/$i --method fbspline --charge-table $record_dir/charge.tab --beta $beta --rcut $cal_rcut --kx $cal_Kx --ky $cal_Ky --kz $cal_Kz --order $cal_order --output-force force.spme.ana.out &>> $mylog

    Lx=`tail $record_dir/$i -n 1 | awk '{printf $1}'`
    Ly=`tail $record_dir/$i -n 1 | awk '{printf $2}'`
    Lz=`tail $record_dir/$i -n 1 | awk '{printf $3}'`
    hhx=`echo "$Lx / $cal_Kx / 2.0" | bc -l`
    hhy=`echo "$Ly / $cal_Ky / 2.0" | bc -l`
    hhz=`echo "$Lz / $cal_Kz / 2.0" | bc -l`
    ./tools/real.error/spme -x $hhx -y $hhy -z $hhz --config $record_dir/$i --charge-table $record_dir/charge.tab --method fbspline --beta $beta --rcut $cal_rcut --kx $cal_Kx --ky $cal_Ky --kz $cal_Kz --order $cal_order --output-force shift.force.spme.ana.out &>> $mylog
    ./tools/real.error/average.force --force-0 force.spme.ana.out --force-1 shift.force.spme.ana.out -o force.spme.ana.st.out
    
    ./tools/real.error/compare.force --config $record_dir/$i --charge-table $record_dir/charge.tab --refh $real_h --force-0 $record_dir/force.ref.$iname.out --force-1 force.spme.ana.st.out -o error.spme.ana.st.$iname.out -m meanf.spme.ana.st.$iname.out &>> $mylog

    mv error.spme.ana.st.$iname.out meanf.spme.ana.st.$iname.out $errors_dir
    echo "$errors_dir/error.spme.ana.st.$iname.out $errors_dir/meanf.spme.ana.st.$iname.out" >> $spme_ana_st_record
done
mv -f $spme_ana_st_record $errors_dir

echo "# #############################################################################" &>> $mylog
echo "# combine results" &>> $mylog
echo "# #############################################################################" &>> $mylog
./tools/real.error/combine.results --record-file $errors_dir/$spme_ana_st_record -o $errors_dir/real.spme.ana.st.error.out -m $errors_dir/real.spme.ana.st.meanf.out &>> $mylog

rm -f force.spme.ana.out shift.force.spme.ana.out force.spme.ana.st.out
mv -f make.log $mylog $errors_dir
