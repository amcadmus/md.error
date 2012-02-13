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
rm -f $errors_dir/parameters.ewald.esti.sh
cp parameters.sh $errors_dir/parameters.ewald.esti.sh

mylog=ewald.esti.log
rm -f $mylog
touch $mylog

make -C ./tools/analyze/ &>> make.log

# esti the rec error
if echo $project_name | grep one_peak &> /dev/null; then
    ./tools/analyze/error.ewald -t $record_dir/traj.xtc -q $record_dir/charge.tab --my-charge $charge --kmaxx $cal_Kx --kmaxy $cal_Ky --kmaxz $cal_Kz --kx $ref_Kx --ky $ref_Ky --kz $ref_Kz --beta $beta &>> $mylog
else
    ./tools/analyze/error.ewald -t $record_dir/traj.xtc --kmaxx $cal_Kx --kmaxy $cal_Ky --kmaxz $cal_Kz --kx $ref_Kx --ky $ref_Ky --kz $ref_Kz --beta $beta &>> $mylog
fi

mv -f rho.x.avg.out $errors_dir/
mv -f error.out $errors_dir/esti.rec.ewald.error.out
mv -f meanf.out $errors_dir/esti.rec.ewald.meanf.out

# combine rec with dir
./tools/analyze/combine.error --dir-meanf $errors_dir/esti.dir.meanf.out --rec-meanf $errors_dir/esti.rec.ewald.meanf.out  --dir-error $errors_dir/esti.dir.error.out --rec-error $errors_dir/esti.rec.ewald.error.out --output-error $errors_dir/esti.ewald.error.out --output-meanf $errors_dir/esti.ewald.meanf.out &>> $mylog


mv -f make.log $mylog $errors_dir
