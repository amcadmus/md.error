set term post enh eps solid font 14 size 8.5cm,6cm
set out 'error.adapt.eps'

set style line 1 lc 1 lw 2 pt 7
set style line 2 lc 3 lw 2 pt 7
set style line 3 lc 2 lw 2 pt 7

set rmargin 8.5

set xrange [0:150]
set x2range [0:150]
set yrange [1e-6:2e-1]
set logscale y
set y2range [0:1]
set mxtics 2
set ytics nomirror
set y2tics .1
set mytics 10
set my2tics 2
set format y '10^{%L}'
set format y2 '%.1f'
set xlabel 'x'
set ylabel 'RMS error of the force'
set y2label 'Density'

pl\
'a.density.x.out' axes x2y2 w l ls 1 not,\
'a.real.x.out'  every 1 axes x1y1 w l ls 3 t '{/Symbol e}_{real}',\
'a.error.x.out' every 1 axes x1y1 w l ls 2 t '{/Symbol e}_{esti}'

