#!/bin/bash

nframe_real_error=50
nframe_traj=100
config_pool=~/study/longRange.hetero/profile.test/config.pool
water_template=$config_pool/water.one_peak

# gen conf parameters
Lx=14.896480
Ly=7.448240
Lz=7.448240
layerWith=0.90
peak_size=7.6
rhoh=32.784
rhol=0.005
charge=0.834

# project command
project_name=uniform
profile_command="$project_name -x $Lx -y $Ly -z $Lz --rho $rhoh -n 0"

peak_dist=5.0
project_name=two_peaks
profile_command="$project_name -d $peak_dist -x $Lx -y $Ly -z $Lz -p $peak_size -t $layerWith -u $rhoh -l $rhol -n 0"

project_name=two_peaks_sep
profile_command="$project_name -x $Lx -y $Ly -z $Lz --layer $layerWith --rhoh $rhoh --rhol $rhol -n 0"

project_name=one_peak
profile_command="$project_name -x $Lx -y $Ly -z $Lz -p $peak_size -t $layerWith -u $rhoh -l $rhol -n 0 --charge $charge"

project_name=water

project_name=rand_water_1
profile_command="$project_name -x $Lx -y $Ly -z $Lz -p $peak_size -t $layerWith -u $rhoh -l $rhol -n 0 --charge $charge"

# param for reference force
beta=2.5
ref_rcut=2.1
ref_Kx=161
ref_Ky=81
ref_Kz=81
# regulate parameters. DO NOT MOVE
nframe_real_error=`printf "%03d" $nframe_real_error`
beta=`printf "%.3f" $beta`
ref_rcut=`printf "%.2f" $ref_rcut`
ref_Kx=`printf "%03d" $ref_Kx`
ref_Ky=`printf "%03d" $ref_Ky`
ref_Kz=`printf "%03d" $ref_Kz`
Lx=`printf "%.1f" $Lx`
Ly=`printf "%.1f" $Ly`
Lz=`printf "%.1f" $Lz`
charge=`printf "%.3f" $charge`
record_dir=$project_name.charge${charge}.box${Lx}x${Ly}x${Lz}.b$beta.r$ref_rcut.K${ref_Kx}x${ref_Ky}x${ref_Kz}.nconf${nframe_real_error}
record_dir=$config_pool/$record_dir

# param for error estimate
cal_rcut=1.3
cal_Kx=121
cal_Ky=61
cal_Kz=61
cal_order=6
real_h=0.5

# regulate parameters. DO NOT MOVE
cal_rcut=`printf "%.2f" $cal_rcut`
cal_Kx=`printf "%03d" $cal_Kx`
cal_Ky=`printf "%03d" $cal_Ky`
cal_Kz=`printf "%03d" $cal_Kz`
cal_order=`printf "%01d" $cal_order`
errors_dir=error.$project_name.charge${charge}.box${Lx}x${Ly}x${Lz}.b$beta.r$cal_rcut.n$cal_order.K${cal_Kx}x${cal_Ky}x${cal_Kz}
