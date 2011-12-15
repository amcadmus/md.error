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
    sed -e "/IndexType densityProfileSamplingFeq/s/=.*/= $densityProfileSamplingFeq;/g"|\
    sed -e "/double refh/s/=.*/= $refh;/g"|\
    sed -e "/IndexType assignForceCorrFeq/s/=.*/= $assignForceCorrFeq;/g"|\
    sed -e "/IndexType updateForceCorrFeq/s/=.*/= $updateForceCorrFeq;/g"|\
    sed -e "/^\#define NThreadsPerBlockCell/s/[0-9][0-9]*/$NTCell/g" |\
    sed -e "/^\#define NThreadsPerBlockAtom/s/[0-9][0-9]*/$NTAtom/g" > tmp
    mv -f tmp $file
}


set_parameter tools.collision/simulation.fcorr/app/lj.nve.cu
make -j4 -C tools.collision/simulation.fcorr/ >& make.log

rm -f conf.gro
if test -f $collsion_init_conf; then
    echo "# the init ball conf exists, generate collision conf."
    make -C tools.collision/gen.conf/ -j4 &>> make.log
    tools.collision/gen.conf/collision.conf -x $collision_x -u $collision_u --ball-conf $collision_init_conf
    mv confout.gro conf.gro
else
    echo "# the init ball conf does not exist, exit"
    exit
fi
natom=`head conf.gro -n 2 | tail -n 1`
natom2=`echo "$natom/2" | bc`
sed -e "s/mola.*/mola $natom2/g" gromacs/topol.top | sed -e "s/molb.*/molb $natom2/g" > tmp.top
mv -f tmp.top gromacs/topol.top
cd gromacs
ln -s ../conf.gro .
ln -s ../traj.xtc .
cd ..
workdir=`pwd`
command="$workdir/tools.collision/simulation.fcorr/lj.nve conf.gro $collision_nstep $device"
echo "command is: $command > gpu-md.out 2> gpu-md.perf"
$command > gpu-md.out 2> gpu-md.perf



