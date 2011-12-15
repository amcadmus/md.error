#!/bin/bash

source parameters.sh

st_start=5000
print_format="%.3f"

step_start=`echo $st_start / $dt | bc -l | cut -d '.' -f 1`
nstart=`echo "$step_start / $thermoFeq" | bc `
echo "# start step is $step_start"

grep -v '#' gpu-md.out | grep -v write > tmp.st
nlines=`wc -l tmp.st | awk '{print $1}'`
ntails=`echo "$nlines - $nstart" | bc `
tail -n $ntails tmp.st > gpu-md.st
rm -f tmp.st

tension_line=`avg_jk  -v col=11 gpu-md.st | grep -v '#'`
tension=`echo $tension_line | awk '{print $1}'`
tension_error=`echo $tension_line | awk '{print $2}'`
tension_error=`echo "$tension_error * 2" | bc -l`

tension_corr=`avg_jk  -v col=14 gpu-md.st | grep -v '#' | awk '{print $1}'`
tension_corr_print=`printf $print_format $tension_corr`

tension_total=`echo "$tension_corr + $tension" | bc -l`
tension_total_print=`printf $print_format $tension_total`

tension_print=`printf $print_format $tension`
tension_error_print=`printf $print_format $tension_error`
echo "# tension error correction corrected"
echo "$tension_print $tension_error_print $tension_corr_print $tension_total_print"


