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
rm -f $errors_dir/parameters.dir.esti.sh
cp parameters.sh $errors_dir/parameters.dir.esti.sh

mylog=dir.esti.log
rm -f $mylog
touch $mylog

make -j8 -C ./tools/analyze/ &>> make.log

# esti the dir error
if echo $project_name | grep one_peak &> /dev/null; then
    ./tools/analyze/error.dir -t $record_dir/traj.xtc -q $record_dir/charge.tab --refh $real_h --beta $beta --rcut $cal_rcut &>> $mylog
else
    ./tools/analyze/error.dir -t $record_dir/traj.xtc --refh $real_h --beta $beta --rcut $cal_rcut &>> $mylog
fi
mv -f rho.x.avg.out $errors_dir/
mv -f error.out $errors_dir/esti.dir.error.out
mv -f meanf.out $errors_dir/esti.dir.meanf.out

mv -f make.log $mylog $errors_dir
