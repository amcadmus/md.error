#!/bin/bash

source parameters.sh

mylog=spme.ana.esti.log
rm -f $mylog
touch $mylog

make -C ./tools/analyze/ &>> make.log

# # esti the dir error
# ./tools/analyze/error.dir -t $record_dir/traj.xtc --refh $real_h --beta $beta --rcut $cal_rcut &>> $mylog
# mv -f rho.x.avg.out $record_dir/
# mv -f error.out $record_dir/esti.dir.error.out
# mv -f meanf.out $record_dir/esti.dir.meanf.out

# esti the rec error
./tools/analyze/error.spme.ana -t $record_dir/traj.xtc -k $cal_K --beta $beta --order $cal_order &>> $mylog
mv -f rho.x.avg.out $record_dir/
mv -f error.out $record_dir/esti.rec.spme.ana.error.out
mv -f meanf.out $record_dir/esti.rec.spme.ana.meanf.out

# combine rec with dir
./tools/analyze/combine.error --dir-meanf $record_dir/esti.dir.meanf.out --rec-meanf $record_dir/esti.rec.spme.ana.meanf.out  --dir-error $record_dir/esti.dir.error.out --rec-error $record_dir/esti.rec.spme.ana.error.out --output-error $record_dir/esti.spme.ana.error.out --output-meanf $record_dir/esti.spme.ana.meanf.out &>> $mylog


mv -f $mylog $record_dir
