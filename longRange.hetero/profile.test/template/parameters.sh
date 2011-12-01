#!/bin/bash

nframe_real_error=5
nframe_error_esti=100

# gen conf parameters
Lx=30
Ly=30
Lz=30
layerWith=2
peak_size=20
rhoh=0.50
rhol=0.05

project_name=one.peak
profile_command="$project_name -x $Lx -y $Ly -z $Lz -p $peak_size -t $layerWith -u $rhoh -l $rhol -n 0"

peak_dist=5.0
project_name=two.peaks
profile_command="$project_name -d $peak_dist -x $Lx -y $Ly -z $Lz -p $peak_size -t $layerWith -u $rhoh -l $rhol -n 0"

project_name=uniform
profile_command="$project_name -x $Lx -y $Ly -z $Lz --rho $rhoh -n 0"

# Ewald parameters
beta=1.0
ref_rcut=5.0
ref_K=65
# regulate parameters. DO NOT MOVE
beta=`printf "%.3f" $beta`
ref_rcut=`printf "%.2f" $ref_rcut`
ref_K=`printf "%03d" $ref_K`
record_dir=$project_name.b$beta.r$ref_rcut.K$ref_K

# naive summation
cal_rcut=3.0
cal_K=59
cal_order=6

real_h=1.0

