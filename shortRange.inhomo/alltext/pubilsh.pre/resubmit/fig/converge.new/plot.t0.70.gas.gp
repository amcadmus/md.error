set term post eps enh color font 16 size 8.5cm,6cm
set out 't0.70.gas.eps'

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
#set yrange [-0.001:0.002]
set xtics 1
set ytics 0.0005
set mxtics 2
set mytics 5
set format y "%.4f"
set xlabel "r_c*"
set ylabel "{/Symbol r}*"

set key right top
set title 'T*=0.70'
unset title
unset label
set label 'T*=0.70, gas' at 5.5, 0.0033
set style fill pattern 8
set clip two

pl \
't0.70.gas.ref.out' u 1:($2+$3):($2-$3) w filledcu lc 0 not,\
't0.70.gas.ref.out' u 1:2 w l lt 1 lc 0 lw 3 not,\
't0.70.gas.uni.out'   u 1:($2):3 w e ls 11 not, '' u 1:($2) w l ls 11 t'URC',\
't0.70.gas.adapt.out' u 1:($2):3 w e ls 21 not, '' u 1:($2) w l ls 21 t'ARC',\
't0.70.gas.fcorr.out' u 1:($2):3 w e ls 31 not, '' u 1:($2) w l ls 31 t'LFC'
