set term post enh eps font 14 size 8.5cm,6cm
set out 'rcut-and-error.eps'

set style line 1  lt 1 lc 1 lw 3 pt 7
set style line 2  lt 1 lc 3 lw 3 pt 7
set style line 21 lt 2 lc 3 lw 3 pt 7
set style line 22 lt 4 lc 3 lw 3 pt 7
set style line 3  lt 1 lc 2 lw 3 pt 5

set rmargin 8.5

set xrange [0:150]
set x2range [0:150]
set yrange [1e-6:2e-1]
set logscale y
set y2range [2:9]
set mxtics 2
set ytics nomirror
set y2tics 1.0
set mytics 10
set my2tics 2
set format y '10^{%L}'
set format y2 '%.1f'
set xlabel 'x'
set ylabel 'RMS error of the force'
set y2label 'Cut-off radius'

pl\
'a.rcut.x.out' every 1 axes x2y2 w l ls 1 t 'r_c',\
'a.real.x.out'  every 1 axes x1y1 w p ls 3 t '{/Symbol e}_{real}',\
'a.error.x.out' every 1 axes x1y1 w l ls 2 t '{/Symbol e}_{esti}',\
'' u 1:3 axes x1y1 w l ls 21 t '{/Symbol e}_{homo}',\
'' u 1:4 axes x1y1 w l ls 22 t '{/Symbol e}_{inhomo}'

