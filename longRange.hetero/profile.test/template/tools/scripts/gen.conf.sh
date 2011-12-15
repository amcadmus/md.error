#!/bin/bash

source parameters.sh

mylog=gen.conf.log

rm -f $mylog
if test -d $record_dir; then
    echo "# existing record dir $record_dir, mv to backup"
    mv -f $record_dir $record_dir.`date +%s`
fi
rm -fr $record_dir
mkdir -p $record_dir

make -j8 -C ./tools/profile/ &> make.log
make -j8 -C ./tools/real.error/ &>> make.log

seed=`date +%s`

for i in `seq 1 $nframe_traj`;
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

mv -f confFilenames.record $record_dir
./tools/profile/conf2traj -r $record_dir/confFilenames.record -o $record_dir/traj.xtc

for i in `seq $(($nframe_real_error+1)) $nframe_traj`;
do
    iname=`printf "%04d" $i`
    rm -f $record_dir/conf.$iname.gro
done

mv -f make.log $mylog $record_dir
rm -f traj.xtc


