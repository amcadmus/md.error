#!/bin/bash

if test ! -d $record_dir; then
    echo "no record dir $record_dir"
    exit
fi

if test ! -f $errors_dir/esti.dir.error.out; then
    ./tools/scripts/dir.esti.sh
fi  

echo "# ewald real"
./tools/scripts/ewald.real.sh
echo "# ewald esti"
./tools/scripts/ewald.esti.sh
