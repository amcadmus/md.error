#!/bin/bash

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
