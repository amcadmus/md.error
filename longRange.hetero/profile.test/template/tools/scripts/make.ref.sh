#!/bin/bash

source parameters.sh

mylog=make.ref.log

rm -f $mylog

make -C ./tools/profile/ &> make.log
make -C ./tools/real.error/ &>> make.log

for i in `ls $record_dir | grep conf | grep gro`;
do
    iname=`echo $i | cut -d '.' -f 2`
    echo "# #############################################################################" &>> $mylog
    echo "# doing $iname" &>> $mylog
    echo "# #############################################################################" &>> $mylog
    ./tools/real.error/ewaldSum --config $record_dir/$i --beta $beta --rcut $ref_rcut --grid-number $ref_K --output-force force.ref.$iname.out &>> $mylog
    mv force.ref.$iname.out $record_dir
done

mv -f $mylog $record_dir
