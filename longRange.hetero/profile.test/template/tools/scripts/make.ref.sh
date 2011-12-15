#!/bin/bash

source parameters.sh

if test ! -d $record_dir; then
    echo "no record dir $record_dir"
    exit
fi
rm -f $record_dir/parameters.make.ref.sh
cp parameters.sh $record_dir/parameters.make.ref.sh

mylog=make.ref.log

rm -f $mylog

make -j8 -C ./tools/profile/ &> make.log
make -j8 -C ./tools/real.error/ &>> make.log

for i in `ls $record_dir | grep conf | grep gro`;
do
    iname=`echo $i | cut -d '.' -f 2`
    echo "# #############################################################################" &>> $mylog
    echo "# doing $iname" &>> $mylog
    echo "# #############################################################################" &>> $mylog
    ./tools/real.error/ewaldSum --config $record_dir/$i --beta $beta --rcut $ref_rcut --kx $ref_Kx --ky $ref_Ky --kz $ref_Kz --output-force force.ref.$iname.out &>> $mylog
    mv force.ref.$iname.out $record_dir
done

mv -f make.log $mylog $record_dir
