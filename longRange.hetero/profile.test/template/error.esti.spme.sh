#!/bin/bash

source parameters.sh

if test ! -d $record_dir; then
    echo "no record dir $record_dir"
    exit
fi

if test ! -f $errors_dir/esti.dir.error.out; then
    ./tools/scripts/dir.esti.sh
fi  

./tools/scripts/spme.ik.real.sh
./tools/scripts/spme.ik.esti.sh
./tools/scripts/spme.ana.real.sh
./tools/scripts/spme.ana.esti.sh

gp_file_names="pl.ana.error.gp  pl.ana.meanf.gp  pl.ik.error.gp  pl.ik.meanf.gp"
for gp_file_name in $gp_file_names; do
#    echo "test ! -f $errors_dir/$gp_file_name && cp tools/gp.scripts/$gp_file_name $errors_dir"
    test ! -f $errors_dir/$gp_file_name && cp tools/gp.scripts/$gp_file_name $errors_dir
done
