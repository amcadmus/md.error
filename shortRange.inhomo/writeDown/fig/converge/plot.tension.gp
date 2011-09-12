set term post eps enh color font 16 size 8.5cm,6cm

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
set mytics 5
set format y "%.2f"
set key bottom right

set out 'tension.t1.20.eps'
set title 'T*=1.20'
set yrange [0.00:0.20]
pl 0.128 w l lc 0 lw 3,\
   't1.20.tension.uni.out' u 1:2:3 w e ls 11 not, '' u 1:2 w l ls 11 t 'URC HPC',\
   '' u 1:($2-$4):3 w e ls 11 not, '' u 1:($2-$4) w l ls 12 t 'URC',\
   't1.20.tension.adapt.out' u 1:2:3 w e ls 21 not, '' u 1:2 w l ls 21 t 'ARC HPC',\
   '' u 1:($2-$4):3 w e ls 21 not, '' u 1:($2-$4) w l ls 22 t 'ARC',\
   't1.20.tension.fcorr.out' u 1:2:3 w e ls 31 not, '' u 1:2 w l ls 31 t 'FC HPC',\
   '' u 1:($2-$4):3 w e ls 31 not, '' u 1:($2-$4) w l ls 32 t 'FC'

set out 'tension.t1.10.eps'
set title 'T*=1.10'
set yrange [0.10:0.40]
pl 0.307 w l lc 0 lw 3,\
   't1.10.tension.uni.out' u 1:2:3 w e ls 11 not, '' u 1:2 w l ls 11 t 'URC HPC',\
   '' u 1:($2-$4):3 w e ls 11 not, '' u 1:($2-$4) w l ls 12 t 'URC',\
   't1.10.tension.adapt.out' u 1:2:3 w e ls 21 not, '' u 1:2 w l ls 21 t 'ARC HPC',\
   '' u 1:($2-$4):3 w e ls 21 not, '' u 1:($2-$4) w l ls 22 t 'ARC',\
   't1.10.tension.fcorr.out' u 1:2:3 w e ls 31 not, '' u 1:2 w l ls 31 t 'FC HPC',\
   '' u 1:($2-$4):3 w e ls 31 not, '' u 1:($2-$4) w l ls 32 t 'FC'

set out 'tension.t0.85.eps'
set title 'T*=0.85'
set yrange [0.50:0.90]
pl 0.828 w l lc 0 lw 3,\
   't0.85.tension.uni.out' u 1:2:3 w e ls 11 not, '' u 1:2 w l ls 11 t 'URC HPC',\
   '' u 1:($2-$4):3 w e ls 11 not, '' u 1:($2-$4) w l ls 12 t 'URC',\
   't0.85.tension.adapt.out' u 1:2:3 w e ls 21 not, '' u 1:2 w l ls 21 t 'ARC HPC',\
   '' u 1:($2-$4):3 w e ls 21 not, '' u 1:($2-$4) w l ls 22 t 'ARC',\
   't0.85.tension.fcorr.out' u 1:2:3 w e ls 31 not, '' u 1:2 w l ls 31 t 'FC HPC',\
   '' u 1:($2-$4):3 w e ls 31 not, '' u 1:($2-$4) w l ls 32 t 'FC'

set out 'tension.t0.70.eps'
set title 'T*=0.70'
set yrange [0.80:1.30]
pl 1.150 w l lc 0 lw 3,\
   't0.70.tension.uni.out' u 1:2:3 w e ls 11 not, '' u 1:2 w l ls 11 t 'URC HPC',\
   '' u 1:($2-$4):3 w e ls 11 not, '' u 1:($2-$4) w l ls 12 t 'URC',\
   't0.70.tension.adapt.out' u 1:2:3 w e ls 21 not, '' u 1:2 w l ls 21 t 'ARC HPC',\
   '' u 1:($2-$4):3 w e ls 21 not, '' u 1:($2-$4) w l ls 22 t 'ARC',\
   't0.70.tension.fcorr.out' u 1:2:3 w e ls 31 not, '' u 1:2 w l ls 31 t 'FC HPC',\
   '' u 1:($2-$4):3 w e ls 31 not, '' u 1:($2-$4) w l ls 32 t 'FC'
   