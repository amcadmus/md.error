#!/bin/bash

source parameters.sh

if test ! -d $record_dir; then
    echo "no record dir $record_dir"
    exit
fi

if echo $project_name | grep -v rand | grep water &> /dev/null; then
    echo "# real spme ik noHO"
    ./tools/scripts/spme.ik.noHO.real.sh
    echo "# real spme ik noHO st"
    ./tools/scripts/spme.ik.noHO.st.real.sh

    echo "# real spme ik"
    ./tools/scripts/spme.ik.real.sh
    echo "# real spme ik st"
    ./tools/scripts/spme.ik.st.real.sh
else
    echo "# real spme ik"
    ./tools/scripts/spme.ik.real.sh
    echo "# real spme ik st"
    ./tools/scripts/spme.ik.st.real.sh
fi

if echo $project_name | grep -v rand | grep water &> /dev/null; then
    echo "# real spme ana noHO"
    ./tools/scripts/spme.ana.noHO.real.sh
    echo "# real spme ana noHO st"
    ./tools/scripts/spme.ana.noHO.st.real.sh

    echo "# real spme ana"
    ./tools/scripts/spme.ana.real.sh
    echo "# real spme ana selfCorr"
    ./tools/scripts/spme.ana.selfCorr.real.sh
    echo "# real spme ana st"
    ./tools/scripts/spme.ana.st.real.sh
else
    echo "# real spme ana"
    ./tools/scripts/spme.ana.real.sh
    echo "# real spme ana selfCorr"
    ./tools/scripts/spme.ana.selfCorr.real.sh
    echo "# real spme ana st"
    ./tools/scripts/spme.ana.st.real.sh
fi

gp_file_names="pl.water.ana.error.gp pl.water.ik.error.gp pl.ana.error.gp  pl.ana.meanf.gp  pl.ik.error.gp  pl.ik.meanf.gp"
for gp_file_name in $gp_file_names; do
#    echo "test ! -f $errors_dir/$gp_file_name && cp tools/gp.scripts/$gp_file_name $errors_dir"
    test ! -f $errors_dir/$gp_file_name && cp tools/gp.scripts/$gp_file_name $errors_dir
done
