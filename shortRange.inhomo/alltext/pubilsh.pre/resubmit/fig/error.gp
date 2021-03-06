set term post eps enh font 14 size 8.5cm,6cm
set out 'error.uniform.eps'

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
't0.85-n16000-rc07.5uni/a.density.x.out' axes x2y2 w l ls 1 t '{/Symbol r}',\
't0.85-n16000-rc07.5uni/a.real.x.out'  every 1 axes x1y1 w p ls 3 t '{/Symbol e}_{real}',\
't0.85-n16000-rc07.5uni/a.error.xtc.x.out' every 1 axes x1y1 w l ls 2 t '{/Symbol e}_{esti}',\
't0.85-n16000-rc07.5uni/a.error.xtc.x.part1.out' every 1 axes x1y1 w l ls 21 t '{/Symbol e}_{homo}',\
't0.85-n16000-rc07.5uni/a.error.xtc.x.part2.out' every 1 axes x1y1 w l ls 22 t '{/Symbol e}_{inhomo}'


set format y '%L'
set ylabel 'RMS error of the force log_{10}[{/Symbol e}*]'
set out 'pre.error.uniform.eps'
set xlabel 'x*'

rep
