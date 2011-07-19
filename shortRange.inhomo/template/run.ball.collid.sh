#!/bin/bash

source p.collision.sh

function set_parameter () {
    file=$1
    sed -e "/IndexType nstep/s/=.*/= $collision_nstep;/g" $file |\
    sed -e "/IndexType confFeq/s/=.*/= $collision_confFeq;/g" |\
    sed -e "/IndexType thermoFeq/s/=.*/= $collision_thermoFeq;/g"|\
    sed -e "/ScalorType rcut/s/=.*/= $collision_rcut;/g" |\
    sed -e "/ScalorType nlistExten/s/=.*/= $collision_nlistExten;/g" |\
    sed -e "/ScalorType nlistSizeFactor/s/=.*/= $collision_nlistSizeFactor;/g" |\
    sed -e "/IndexType  clistDivision/s/=.*/= $collision_clistDivision;/g" |\
    sed -e "/ScalorType dt/s/=.*/= $collision_dt;/g" |\
    sed -e "/^\#define NThreadsPerBlockCell/s/[0-9][0-9]*/$NTCell/g" |\
    sed -e "/^\#define NThreadsPerBlockAtom/s/[0-9][0-9]*/$NTAtom/g" > tmp
    mv -f tmp $file
}


set_parameter tools/simulation.uniform.cutoff/app/lj.nve.cu
make -j4 -C tools/simulation.uniform.cutoff/ &> make.log

rm -f conf.gro
if test -f $collsion_init_conf; then
    echo "# the init ball conf exists, generate collision conf."
    make -C tools.collision/gen.conf/ -j4 &>> make.log
    tools.collision/gen.conf/collision.conf -x $collision_x -u $collision_u
    mv confout.gro conf.gro
else
    echo "# the init ball conf does not exist, exit"
    exit
fi
workdir=`pwd`
command="$workdir/tools/simulation.uniform.cutoff/lj.nve conf.gro $collision_nstep $device"
echo "command is: $command > gpu-md.out 2> gpu-md.perf"
$command > gpu-md.out 2> gpu-md.perf



