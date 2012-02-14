#!/bin/bash

nframe_real_error=1
nframe_traj=100
config_pool=~/study/longRange.hetero/profile.test/config.pool

# gen conf parameters
Lx=40
Ly=20
Lz=20
layerWith=2
peak_size=20
rhoh=1.00
rhol=0.05
charge=1.0

project_name=uniform
profile_command="$project_name -x $Lx -y $Ly -z $Lz --rho $rhoh -n 0"

peak_dist=5.0
project_name=two_peaks
profile_command="$project_name -d $peak_dist -x $Lx -y $Ly -z $Lz -p $peak_size -t $layerWith -u $rhoh -l $rhol -n 0"

project_name=two_peaks_sep
profile_command="$project_name -x $Lx -y $Ly -z $Lz --layer $layerWith --rhoh $rhoh --rhol $rhol -n 0"

project_name=one_peak
profile_command="$project_name -x $Lx -y $Ly -z $Lz -p $peak_size -t $layerWith -u $rhoh -l $rhol -n 0 --charge $charge"

# param for reference force
beta=1.0
ref_rcut=6.0
ref_Kx=129
ref_Ky=65
ref_Kz=65
# regulate parameters. DO NOT MOVE
beta=`printf "%.3f" $beta`
ref_rcut=`printf "%.2f" $ref_rcut`
ref_Kx=`printf "%03d" $ref_Kx`
ref_Ky=`printf "%03d" $ref_Ky`
ref_Kz=`printf "%03d" $ref_Kz`
charge=`printf "%.3f" $charge`
if echo $project_name | grep one_peak &> /dev/null; then
    record_dir=$project_name.charge${charge}.box${Lx}x${Ly}x${Lz}.b$beta.r$ref_rcut.K${ref_Kx}x${ref_Ky}x${ref_Kz}
else
    record_dir=$project_name.box${Lx}x${Ly}x${Lz}.b$beta.r$ref_rcut.K${ref_Kx}x${ref_Ky}x${ref_Kz}
fi
record_dir=$config_pool/$record_dir

# param for error estimate
cal_rcut=3.0
cal_Kx=101
cal_Ky=51
cal_Kz=51
cal_order=6
real_h=1.0

# regulate parameters. DO NOT MOVE
cal_rcut=`printf "%.2f" $cal_rcut`
cal_Kx=`printf "%03d" $cal_Kx`
cal_Ky=`printf "%03d" $cal_Ky`
cal_Kz=`printf "%03d" $cal_Kz`
cal_order=`printf "%01d" $cal_order`
if echo $project_name | grep one_peak &> /dev/null; then
    errors_dir=error.$project_name.charge${charge}.box${Lx}x${Ly}x${Lz}.b$beta.r$cal_rcut.n$cal_order.K${cal_Kx}x${cal_Ky}x${cal_Kz}
else
    errors_dir=error.$project_name.box${Lx}x${Ly}x${Lz}.b$beta.r$cal_rcut.n$cal_order.K${cal_Kx}x${cal_Ky}x${cal_Kz}
fi
