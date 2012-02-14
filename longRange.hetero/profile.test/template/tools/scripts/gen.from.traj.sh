#!/bin/bash

source parameters.sh

mylog=gen.from.traj.log

rm -f $mylog
if test ! -f $record_dir/traj.xtc; then
    echo "# no file $record_dir/traj.xtc"
    exit
fi
if test ! -f $record_dir/charge.tab; then
    echo "# no file $record_dir/charge.tab"
    exit
fi
cp parameters.sh $record_dir/parameters.gen.from.traj.sh

make -j8 -C ./tools/profile/ &> make.log

./tools/profile/split.traj -i $record_dir/traj.xtc &>> $mylog
mv -f conf*gro $record_dir/

cd $record_dir/
rm -fr tmp
mkdir tmp
mv `ls | grep conf | grep gro | head -n $nframe_real_error ` tmp
rm -f conf*gro
mv tmp/*gro .
rm -fr tmp
cd -

mv -f make.log $mylog $record_dir


