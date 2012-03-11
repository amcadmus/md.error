#!/bin/bash

source parameters.sh

if test ! -d $water_template; then
    echo "no water template $water_template"
    exit
fi

if echo $project_name | grep water | grep -v rand &> /dev/null; then
    if test -d $record_dir; then
	echo "# find $record_dir, backup"
	mv $record_dir $record_dir.`date +%s`
    fi
    cp -a $water_template $record_dir
    echo "# gen conf from traj"
    ./tools/scripts/gen.from.traj.sh
    echo "# make ref"
#    ./tools/scripts/make.ref.sh
    ./tools/scripts/make.noHO.ref.sh
else
    echo "# gen conf"
    ./tools/scripts/gen.conf.sh
    echo "# make ref"
    ./tools/scripts/make.ref.sh
fi

#echo "# make ref"
#./tools/scripts/make.ref.sh
#./tools/scripts/make.noHO.ref.sh


