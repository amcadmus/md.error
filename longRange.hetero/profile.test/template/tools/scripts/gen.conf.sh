#!/bin/bash

source parameters.sh

mylog=gen.conf.log

rm -f $mylog
rm -fr $record_dir
mkdir -p $record_dir

make -C ./tools/profile/ &> make.log
make -C ./tools/real.error/ &>> make.log

for i in `seq 1 $nframe_real_error`;
do
    iname=`printf "%04d" $i`
    echo "# #############################################################################" &>> $mylog
    echo "# doing $iname" &>> $mylog
    echo "# #############################################################################" &>> $mylog
    ./tools/profile/$profile_command -s `date +%s` &>> $mylog
    mv -f conf.gro $record_dir/conf.$iname.gro
    echo "conf.$iname.gro" >> confFilenames.record
done

mv -f $mylog confFilenames.record $record_dir

