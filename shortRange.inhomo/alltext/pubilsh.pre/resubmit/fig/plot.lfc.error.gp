set term post enh eps font 14 size 8.5cm,6cm
set out 'fcorr.and.error.eps'

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
# set y2range [-0.005:0.005]
set mxtics 2
set ytics nomirror
set y2tics 1
set mytics 10
set my2tics 2
set format y '10^{%L}'
set format y2 '%.f'
set xlabel 'x'
set ylabel 'RMS error of the force'
set y2label 'Mean error force ( x 10^3)'


pl\
't0.85-n16000-fcorr-rc07.5-feq0200/afc.x.out' u 1:($2*1e3) every 1 axes x1y2 w l ls 1 t '<{/Symbol D}F_{x}>',\
't0.85-n16000-fcorr-rc07.5-feq0200/a.real.x.out'  every 1 axes x1y1 w p ls 3 t '{/Symbol e}_{real}',\
't0.85-n16000-fcorr-rc07.5-feq0200/a.error.xtc.x.part1.out' axes x1y1 w l ls 2 t '{/Symbol e}_{homo}'

set format y '%L'
set xlabel 'x*'
set ylabel 'RMS force error  log_{10} ({/Symbol e}*)'
set y2label 'Mean error force  10^3<{/Symbol D}F_x*>'
set out 'pre.fcorr.and.error.eps'

rep

