#!/bin/bash

source parameters.sh

mylog=ewald.esti.log
rm -f $mylog
touch $mylog

make -C ./tools/analyze/ &>> make.log

# # esti the dir error
# ./tools/analyze/error.dir -t $record_dir/traj.xtc --refh $real_h --beta $beta --rcut $cal_rcut &>> $mylog
# mv -f rho.x.avg.out $record_dir/
# mv -f error.out $record_dir/esti.dir.error.out
# mv -f meanf.out $record_dir/esti.dir.meanf.out

# esti the rec error
./tools/analyze/error.ewald -t $record_dir/traj.xtc --kmax $cal_K -k $ref_K --beta $beta &>> $mylog
mv -f rho.x.avg.out $record_dir/
mv -f error.out $record_dir/esti.rec.ewald.error.out
mv -f meanf.out $record_dir/esti.rec.ewald.meanf.out

mv -f $mylog $record_dir
