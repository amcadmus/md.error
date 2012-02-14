#!/bin/bash

source parameters.sh

if test ! -d $record_dir; then
    echo "no record dir $record_dir"
    exit
fi

if test ! -f $errors_dir/esti.dir.error.out; then
    echo "# esti spme dir"
    ./tools/scripts/dir.esti.sh
fi  

echo "# real spme ik"
./tools/scripts/spme.ik.real.sh
echo "# esti spme ik"
./tools/scripts/spme.ik.esti.sh
echo "# real spme ana"
./tools/scripts/spme.ana.real.sh
echo "# esti spme ana"
./tools/scripts/spme.ana.esti.sh

gp_file_names="pl.ana.error.gp  pl.ana.meanf.gp  pl.ik.error.gp  pl.ik.meanf.gp"
for gp_file_name in $gp_file_names; do
#    echo "test ! -f $errors_dir/$gp_file_name && cp tools/gp.scripts/$gp_file_name $errors_dir"
    test ! -f $errors_dir/$gp_file_name && cp tools/gp.scripts/$gp_file_name $errors_dir
done
