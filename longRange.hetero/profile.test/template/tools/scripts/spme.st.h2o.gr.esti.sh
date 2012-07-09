#!/bin/bash

source parameters.sh

if test ! -d $record_dir; then
    echo "no record dir $record_dir"
    exit
fi
if test ! -d $errors_dir; then
    echo "# no errors dir $errors_dir, make it"
    mkdir -p $errors_dir
fi
rm -f $errors_dir/parameters.spme.st.h2o.gr.esti.sh
cp parameters.sh $errors_dir/parameters.spme.st.h2o.gr.esti.sh

mylog=spme.st.h2o.gr.esti.log
rm -f $mylog
touch $mylog

make -C -j8 ./tools/analyze/ &>> make.log

# esti the rec error
if test ! -f $record_dir/water.gro; then
    echo "# no orientation file $record_dir/water.gro"
    exit
fi
./tools/analyze/error.spme.st.h2o.gr -t $record_dir/traj.xtc -q $record_dir/charge.tab --my-charge $charge --kx $cal_Kx --ky $cal_Ky --kz $cal_Kz --beta $beta --order $cal_order --orientation $record_dir/water.gro --rhoMol $rhoh --gr-up $gr_up --rdf-OO $record_dir/rdf.oo.1e-2.xvg --rdf-OH $record_dir/rdf.oh.1e-2.xvg --rdf-HH $record_dir/rdf.hh.1e-2.xvg &>> $mylog
mv -f rho.x.avg.out $errors_dir/
mv -f error.out $errors_dir/esti.rec.spme.st.h2o.gr.error.out
mv -f meanf.out $errors_dir/esti.rec.spme.st.h2o.gr.meanf.out

# combine rec with dir
./tools/analyze/combine.error --dir-meanf $errors_dir/esti.dir.meanf.out --rec-meanf $errors_dir/esti.rec.spme.st.h2o.gr.meanf.out  --dir-error $errors_dir/esti.dir.error.out --rec-error $errors_dir/esti.rec.spme.st.h2o.gr.error.out --output-error $errors_dir/esti.spme.st.h2o.gr.error.out --output-meanf $errors_dir/esti.spme.st.h2o.gr.meanf.out &>> $mylog

mv -f make.log $mylog $errors_dir
