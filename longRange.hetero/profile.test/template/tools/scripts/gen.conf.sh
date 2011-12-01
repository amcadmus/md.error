#!/bin/bash

source parameters.sh

mylog=gen.conf.log

rm -f $mylog
rm -fr $record_dir
mkdir -p $record_dir

make -j8 -C ./tools/profile/ &> make.log
make -j8 -C ./tools/real.error/ &>> make.log

seed=`date +%s`

for i in `seq 1 $nframe_real_error`;
do
    iname=`printf "%04d" $i`
    echo "# #############################################################################" &>> $mylog
    echo "# doing $iname" &>> $mylog
    echo "# #############################################################################" &>> $mylog
    ./tools/profile/$profile_command -s $seed &>> $mylog
    seed=$(($seed+1))
    mv -f conf.gro $record_dir/conf.$iname.gro
    echo "$record_dir/conf.$iname.gro" >> confFilenames.record
done

mv -f $mylog confFilenames.record $record_dir

./tools/profile/conf2traj -r $record_dir/confFilenames.record -o $record_dir/traj.xtc

