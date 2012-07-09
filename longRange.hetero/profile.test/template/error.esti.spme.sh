#!/bin/bash

source parameters.sh

if test ! -d $record_dir; then
    echo "no record dir $record_dir"
    exit
fi

echo "# esti spme dir"
./tools/scripts/dir.esti.sh
echo "# esti spme ik"
./tools/scripts/spme.ik.esti.sh
echo "# esti spme ana"
./tools/scripts/spme.ana.esti.sh
echo "# esti spme ana self"
./tools/scripts/spme.ana.selfCorr.esti.sh
echo "# esti spme st"
./tools/scripts/spme.st.esti.sh
if echo $project_name | grep -v rand | grep water &> /dev/null; then
    echo "# esti spme ik h2o"
    ./tools/scripts/spme.ik.h2o.esti.sh
    echo "# esti spme ana h2o"
    ./tools/scripts/spme.ana.h2o.esti.sh
    echo "# esti spme st h2o"
    ./tools/scripts/spme.st.h2o.esti.sh
    echo "# esti spme st h2o"
    ./tools/scripts/spme.st.h2o.gr.esti.sh
fi

gp_file_names="pl.water.ana.error.gp pl.water.ik.error.gp pl.ana.error.gp  pl.ana.meanf.gp  pl.ik.error.gp  pl.ik.meanf.gp"
for gp_file_name in $gp_file_names; do
#    echo "test ! -f $errors_dir/$gp_file_name && cp tools/gp.scripts/$gp_file_name $errors_dir"
    test ! -f $errors_dir/$gp_file_name && cp tools/gp.scripts/$gp_file_name $errors_dir
done
