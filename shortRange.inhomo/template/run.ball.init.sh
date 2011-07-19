#!/bin/bash

source p.collision.sh

function set_parameter () {
    file=$1
    sed -e "/IndexType nstep/s/=.*/= $ball_nstep;/g" $file |\
    sed -e "/IndexType confFeq/s/=.*/= $ball_confFeq;/g" |\
    sed -e "/IndexType thermoFeq/s/=.*/= $ball_thermoFeq;/g"|\
    sed -e "/ScalorType rcut/s/=.*/= $ball_rcut;/g" |\
    sed -e "/ScalorType nlistExten/s/=.*/= $ball_nlistExten;/g" |\
    sed -e "/ScalorType nlistSizeFactor/s/=.*/= $ball_nlistSizeFactor;/g" |\
    sed -e "/ScalorType refT/s/=.*/= $ball_refT;/g" |\
    sed -e "/ScalorType tauT/s/=.*/= $ball_tauT;/g" |\
    sed -e "/ScalorType dt/s/=.*/= $ball_dt;/g" |\
    sed -e "/^\#define NThreadsPerBlockCell/s/[0-9][0-9]*/$NTCell/g" |\
    sed -e "/^\#define NThreadsPerBlockAtom/s/[0-9][0-9]*/$NTAtom/g" > tmp
    mv -f tmp $file
}


set_parameter tools/simulation.uniform.cutoff/app/lj.nhc.cu
make -j4 -C tools/simulation.uniform.cutoff/ &> make.log

rm -f conf.gro
if test -f $ball_init_conf; then
    echo "# the init conf exists"
    cp $ball_init_conf ./conf.gro
else
    echo "# the init conf does not exist, generate one"
    make -j4 -C tools.collision/gen.conf/ &> make.log
    tools.collision/gen.conf/fcc.ball
    mv confout.gro conf.gro
fi
workdir=`pwd`
command="$workdir/tools/simulation.uniform.cutoff/lj.nhc conf.gro $ball_nstep $device"
echo "command is: $command > gpu-md.out 2> gpu-md.perf"
$command > gpu-md.out 2> gpu-md.perf
cp confout.gro ball.init.gro



