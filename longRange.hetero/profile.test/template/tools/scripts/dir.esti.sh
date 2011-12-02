#!/bin/bash

source parameters.sh

if test ! -d $record_dir; then
    echo "no record dir $record_dir"
    exit
fi

mylog=dir.esti.log
rm -f $mylog
touch $mylog

make -j8 -C ./tools/analyze/ &>> make.log

# esti the dir error
./tools/analyze/error.dir -t $record_dir/traj.xtc --refh $real_h --beta $beta --rcut $cal_rcut &>> $mylog
mv -f rho.x.avg.out $record_dir/
mv -f error.out $record_dir/esti.dir.error.out
mv -f meanf.out $record_dir/esti.dir.meanf.out

mv -f make.log $mylog $record_dir
