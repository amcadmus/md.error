set term post eps enh color font 16 size 8.5cm,6cm
set out 't1.20.eps'

set style line 11 lt 1 lc 1 lw 3 pt 7
set style line 12 lt 2 lc 1 lw 3 pt 7
set style line 13 lt 4 lc 1 lw 3 pt 7
set style line 21 lt 1 lc 2 lw 3 pt 7
set style line 22 lt 2 lc 2 lw 3 pt 7
set style line 23 lt 4 lc 2 lw 3 pt 7
set style line 31 lt 1 lc 3 lw 3 pt 7
set style line 32 lt 2 lc 3 lw 3 pt 7
set style line 33 lt 4 lc 3 lw 3 pt 7

set xrange [2:12]
set xtics 1
set mxtics 2
set yrange [-0.030:0.030]
set ytics 0.010
set mytics 10
set format y "%.3f"
set xlabel "r_c*"
set ylabel "{/Symbol r}* - {/Symbol r}*_{ref}"

set key right bottom
set title 'T*=1.20'
unset title
unset label
set label 'T*=1.20' at 9, 0.022

pl 0 lc 0 lw 3 lt 2 not,\
't1.20.gas.uni.out'    u 1:($2-$4):3 w e ls 11 not, '' u 1:($2-$4) w l ls 12 t'gas URC',\
't1.20.gas.adapt.out'  u 1:($2-$4):3 w e ls 21 not, '' u 1:($2-$4) w l ls 22 t'gas ARC',\
't1.20.gas.fcorr.out'   u 1:($2-$4):3 w e ls 31 not, '' u 1:($2-$4) w l ls 32 t'gas LFC',\
't1.20.liquid.uni.out'   u 1:($2-$4):3 w e ls 11 not, '' u 1:($2-$4) w l ls 11 t'liquid URC',\
't1.20.liquid.adapt.out' u 1:($2-$4):3 w e ls 21 not, '' u 1:($2-$4) w l ls 21 t'liquid ARC',\
't1.20.liquid.fcorr.out'  u 1:($2-$4):3 w e ls 31 not, '' u 1:($2-$4) w l ls 31 t'liquid LFC'
