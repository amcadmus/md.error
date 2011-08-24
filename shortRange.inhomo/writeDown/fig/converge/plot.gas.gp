set term post eps enh color
set out 'gas.eps'

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

pl 0 lc 0 lw 3 lt 2,\
't0.70.gas.uni.out'   u 1:($2-$4):3 w e ls 11, '' u 1:($2-$4) w l ls 11 not,\
't0.70.gas.adapt.out' u 1:($2-$4):3 w e ls 11, '' u 1:($2-$4) w l ls 12 not,\
't0.70.gas.fcorr.out'  u 1:($2-$4):3 w e ls 11, '' u 1:($2-$4) w l ls 13 not,\
't0.85.gas.uni.out'   u 1:($2-$4):3 w e ls 21, '' u 1:($2-$4) w l ls 21 not,\
't0.85.gas.adapt.out' u 1:($2-$4):3 w e ls 21, '' u 1:($2-$4) w l ls 22 not,\
't0.85.gas.fcorr.out'  u 1:($2-$4):3 w e ls 21, '' u 1:($2-$4) w l ls 23 not,\
't1.10.gas.uni.out'   u 1:($2-$4):3 w e ls 31, '' u 1:($2-$4) w l ls 31 not,\
't1.10.gas.adapt.out' u 1:($2-$4):3 w e ls 31, '' u 1:($2-$4) w l ls 32 not,\
't1.10.gas.fcorr.out'  u 1:($2-$4):3 w e ls 31, '' u 1:($2-$4) w l ls 33 not
