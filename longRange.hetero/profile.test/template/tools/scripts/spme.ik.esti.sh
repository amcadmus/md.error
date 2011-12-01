#!/bin/bash

source parameters.sh

mylog=spme.ik.esti.log
rm -f $mylog
touch $mylog

make -C ./tools/analyze/ &>> make.log

# # esti the dir error
# ./tools/analyze/error.dir -t $record_dir/traj.xtc --refh $real_h --beta $beta --rcut $cal_rcut &>> $mylog
# mv -f rho.x.avg.out $record_dir/
# mv -f error.out $record_dir/esti.dir.error.out
# mv -f meanf.out $record_dir/esti.dir.meanf.out

# esti the rec error
./tools/analyze/error.spme.ik -t $record_dir/traj.xtc -k $cal_K --beta $beta --order $cal_order &>> $mylog
mv -f rho.x.avg.out $record_dir/
mv -f error.out $record_dir/esti.rec.spme.ik.error.out
mv -f meanf.out $record_dir/esti.rec.spme.ik.meanf.out

# combine rec with dir
./tools/analyze/combine.error --dir-meanf $record_dir/esti.dir.meanf.out --rec-meanf $record_dir/esti.rec.spme.ik.meanf.out  --dir-error $record_dir/esti.dir.error.out --rec-error $record_dir/esti.rec.spme.ik.error.out --output-error $record_dir/esti.spme.ik.error.out --output-meanf $record_dir/esti.spme.ik.meanf.out &>> $mylog


mv -f $mylog $record_dir
