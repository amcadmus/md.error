#!/bin/bash

source parameters.sh

if echo $project_prefix | grep water &> /dev/null; then
    if test -d $record_dir; then
	echo "# find $record_dir, backup"
	mv $record_dir $record_dir.`date +%s`
    fi
    cp -a $water_template $record_dir
    echo "# gen conf from traj"
    ./tools/scripts/gen.from.traj.sh
else
    echo "# gen conf"
    ./tools/scripts/gen.conf.sh
fi

echo "# make ref"
./tools/scripts/make.ref.sh


