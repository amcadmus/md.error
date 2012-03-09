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
rm -f $errors_dir/parameters.spme.ana.selfCorr.real.sh
cp parameters.sh $errors_dir/parameters.spme.ana.selfCorr.real.sh

mylog=spme.ana.selfCorr.real.log
rm -f $mylog
spme_ana_selfCorr_record=spme.ana.selfCorr.record
rm -f $spme_ana_selfCorr_record

make -j8 -C ./tools/real.error/ &>> make.log

for i in `ls $record_dir | grep conf | grep gro`;
do
    iname=`echo $i | cut -d '.' -f 2`
    echo "# #############################################################################" &>> $mylog
    echo "# doing $iname" &>> $mylog
    echo "# #############################################################################" &>> $mylog
    echo "command is: ./tools/real.error/spme --config $record_dir/$i --method fbspline --self-corr yes --charge-table $record_dir/charge.tab --beta $beta --rcut $cal_rcut --kx $cal_Kx --ky $cal_Ky --kz $cal_Kz --order $cal_order --output-force force.spme.ana.selfCorr.out"
    ./tools/real.error/spme --config $record_dir/$i --method fbspline --self-corr yes --charge-table $record_dir/charge.tab --beta $beta --rcut $cal_rcut --kx $cal_Kx --ky $cal_Ky --kz $cal_Kz --order $cal_order --output-force force.spme.ana.selfCorr.out &>> $mylog
    tools/real.error/compare.force --config $record_dir/$i --charge-table $record_dir/charge.tab --refh $real_h --force-0 $record_dir/force.ref.$iname.out --force-1 force.spme.ana.selfCorr.out -o error.spme.ana.selfCorr.$iname.out -m meanf.spme.ana.selfCorr.$iname.out &>> $mylog

    mv error.spme.ana.selfCorr.$iname.out meanf.spme.ana.selfCorr.$iname.out $errors_dir
    echo "$errors_dir/error.spme.ana.selfCorr.$iname.out $errors_dir/meanf.spme.ana.selfCorr.$iname.out" >> $spme_ana_selfCorr_record
done
mv -f $spme_ana_selfCorr_record $errors_dir

echo "# #############################################################################" &>> $mylog
echo "# combine results" &>> $mylog
echo "# #############################################################################" &>> $mylog
./tools/real.error/combine.results --record-file $errors_dir/$spme_ana_selfCorr_record -o $errors_dir/real.spme.ana.selfCorr.error.out -m $errors_dir/real.spme.ana.selfCorr.meanf.out &>> $mylog

rm -f force.spme.ana.selfCorr.out
mv -f make.log $mylog $errors_dir
