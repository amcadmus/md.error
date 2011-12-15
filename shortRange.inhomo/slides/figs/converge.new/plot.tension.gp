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

set xlabel "r_c*"
set ylabel "{/Symbol g}*"
set xrange [2:12]
set xtics 1
set mxtics 2
set mytics 5
set format y "%.2f"
set key bottom right
set clip two
set style fill pattern 8


set out 'tension-t0p85.eps'
# set title 'T*=0.85'
unset label
set label 'T*=0.85' at 6.0, 0.95
set yrange [0.50:1.00]
pl \
   't0.85.tension.ref.out' u 1:($2+$3):($2-$3) w filledcu lc 0 not,\
   't0.85.tension.ref.out' u 1:2 w l lt 1 lc 0 lw 3 not,\
   't0.85.tension.uni.out'   u 1:($2-$4):3 w e ls 11 not, '' u 1:($2-$4) w l ls 12 t 'URC',\
   't0.85.tension.adapt.out' u 1:($2-$4):3 w e ls 21 not, '' u 1:($2-$4) w l ls 22 t 'ARC',\
   't0.85.tension.fcorr.out' u 1:($2-$4):3 w e ls 31 not, '' u 1:($2-$4) w l ls 32 t 'LFC',\
   't0.85.tension.uni.out'   u 1:2:3 w e ls 11 not, '' u 1:2 w l ls 11 t 'URC+LPC',\
   't0.85.tension.adapt.out' u 1:2:3 w e ls 21 not, '' u 1:2 w l ls 21 t 'ARC+LPC',\
   't0.85.tension.fcorr.out' u 1:2:3 w e ls 31 not, '' u 1:2 w l ls 31 t 'LFC+LPC'


set out 'tension-t0p85-1.eps'

pl \
   't0.85.tension.ref.out' u 1:($2+$3):($2-$3) w filledcu lc 0 not,\
   't0.85.tension.ref.out' u 1:2 w l lt 1 lc 0 lw 3 not,\
   't0.85.tension.uni.out'   u 1:($2-$4):3 w e ls 11 not, '' u 1:($2-$4) w l ls 11 not
