set term post enh eps solid   font 16 size 9.5cm,6cm
set out 'fcorr.eps'

set style line 1 lc 1 lw 3 pt 7
set style line 2 lc 3 lw 3 pt 7
set style line 3 lc 3 lw 3 pt 7

set xrange [0:150]
set x2range [0:150]
# set yrange [0:10]
set y2range [0:1]
set mxtics 2
set ytics nomirror
set y2tics
set mytics 2
set my2tics 4
set xlabel 'x'
set ylabel 'Correction Force'
set y2label 'Density'

pl\
'a.density.x.out' axes x2y2 w l ls 1 not,\
'afc.x.out' every 1 axes x1y1 w l ls 2 t 'F_{corr,x}'
