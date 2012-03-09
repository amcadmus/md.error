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
rm -f $errors_dir/parameters.spme.ana.esti.sh
cp parameters.sh $errors_dir/parameters.spme.ana.esti.sh

mylog=spme.ana.esti.log
rm -f $mylog
touch $mylog

make -j8 -C ./tools/analyze/ &>> make.log

# esti the rec error
./tools/analyze/error.spme.ana -t $record_dir/traj.xtc -q $record_dir/charge.tab --my-charge $charge --kx $cal_Kx --ky $cal_Ky --kz $cal_Kz --beta $beta --order $cal_order &>> $mylog
mv -f rho.x.avg.out $errors_dir/
mv -f error.out $errors_dir/esti.rec.spme.ana.error.out
mv -f meanf.out $errors_dir/esti.rec.spme.ana.meanf.out

# combine rec with dir
./tools/analyze/combine.error --dir-meanf $errors_dir/esti.dir.meanf.out --rec-meanf $errors_dir/esti.rec.spme.ana.meanf.out  --dir-error $errors_dir/esti.dir.error.out --rec-error $errors_dir/esti.rec.spme.ana.error.out --output-error $errors_dir/esti.spme.ana.error.out --output-meanf $errors_dir/esti.spme.ana.meanf.out &>> $mylog

mv -f make.log $mylog $errors_dir
