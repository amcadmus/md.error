#!/bin/bash

source parameters.sh

st_start=5000
print_format="%.1e"

step_start=`echo $st_start / $dt | bc -l | cut -d '.' -f 1`
nstart=`echo "$step_start / $thermoFeq" | bc `
echo "# start step is $step_start"

grep -v '#' gpu-md.out | grep -v write > tmp.st
nlines=`wc -l tmp.st | awk '{print $1}'`
ntails=`echo "$nlines - $nstart" | bc `
tail -n $ntails tmp.st > gpu-md.st
rm -f tmp.st

cost_line=`avg_jk  -v col=15 gpu-md.st | grep -v '#'`
cost=`echo $cost_line | awk '{print $1}'`
cost_error=`echo $cost_line | awk '{print $2}'`

cost_print=`printf $print_format $cost`
cost_error_print=`printf $print_format $cost_error`

echo "# cost error "
echo "$cost_print $cost_error_print "


