#!/bin/bash

source parameters.sh

if test ! -d $record_dir; then
    echo "no record dir $record_dir"
    exit
fi
rm -f $record_dir/parameters.make.noHO.ref.sh
cp parameters.sh $record_dir/parameters.make.noHO.ref.sh

mylog=make.noHO.ref.log

rm -f $mylog

make -j8 -C ./tools/profile/ &> make.log
make -j8 -C ./tools/real.error/ &>> make.log
make -j8 -C ./tools/real.error.exclude.H2O/ &>> make.log

for i in `ls $record_dir | grep conf | grep gro`;
do
    iname=`echo $i | cut -d '.' -f 2`
    echo "# #############################################################################" &>> $mylog
    echo "# doing $iname" &>> $mylog
    echo "# #############################################################################" &>> $mylog
    ./tools/real.error.exclude.H2O/ewaldSum --config $record_dir/$i --charge-table $record_dir/charge.tab --beta $beta --rcut $ref_rcut --kx $ref_Kx --ky $ref_Ky --kz $ref_Kz --output-force force.noHO.ref.$iname.out &>> $mylog
    mv force.noHO.ref.$iname.out $record_dir
done

mv -f make.log $mylog $record_dir
