#!/bin/bash

source parameters.sh

function set_parameter () {
    file=$1
    sed -e "/IndexType nstep/s/=.*/= $nstep;/g" $file |\
    sed -e "/IndexType confFeq/s/=.*/= $confFeq;/g" |\
    sed -e "/IndexType thermoFeq/s/=.*/= $thermoFeq;/g"|\
    sed -e "/ScalorType rcut/s/=.*/= $rcut;/g" |\
    sed -e "/ScalorType nlistExten/s/=.*/= $nlistExten;/g" |\
    sed -e "/ScalorType refT/s/=.*/= $refT;/g" |\
    sed -e "/ScalorType tauT/s/=.*/= $tauT;/g" |\
    sed -e "/ScalorType dt/s/=.*/= $dt;/g" |\
    sed -e "/IndexType rcutUpdateFeq/s/=.*/= $rcutUpdateFeq;/g"|\
    sed -e "/IndexType rcutAssignFeq/s/=.*/= $rcutAssignFeq;/g"|\
    sed -e "/IndexType rcutNumRefine/s/=.*/= $rcutNumRefine;/g"|\
    sed -e "/IndexType densityProfileSamplingFeq/s/=.*/= $densityProfileSamplingFeq;/g"|\
    sed -e "/double refh/s/=.*/= $refh;/g"|\
    sed -e "/double rcmin/s/=.*/= $rcmin;/g"|\
    sed -e "/double rcmax/s/=.*/= $rcmax;/g"|\
    sed -e "/double rcstep/s/=.*/= $rcstep;/g"|\
    sed -e "/double targetPrec/s/=.*/= $targetPrec;/g"|\
    sed -e "/^\#define NThreadsPerBlockCell/s/[0-9][0-9]*/$NTCell/g" |\
    sed -e "/^\#define NThreadsPerBlockAtom/s/[0-9][0-9]*/$NTAtom/g" > tmp
    mv -f tmp $file
}


set_parameter tools/simulation.hetero.cutoff/app/lj.nhc-xtc.cu
make -j4 -C tools/simulation.hetero.cutoff/ &> make.log
rm -f conf.gro
cp $confFile conf.gro
workdir=`pwd`
command="$workdir/tools/simulation.hetero.cutoff/lj.nhc-xtc conf.gro $nstep $device"
echo "command is: $command > gpu-md.out 2> gpu-md.perf"
$command > gpu-md.out 2> gpu-md.perf



