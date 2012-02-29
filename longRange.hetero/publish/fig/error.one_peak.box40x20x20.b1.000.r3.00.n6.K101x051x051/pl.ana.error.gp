set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig.ana.error.eps'

set style line 11 lc 1 lw 3 pt 7
set style line 21 lc 2 lw 3 pt 7
set style line 31 lc 3 lw 3 pt 7

set format y '%.1e'
set xlabel 'x'
set ylabel 'RMS error'
set yrange [0:.8e-3]
set mxtics 5
set mytics 5

pl\
'real.spme.ana.error.out' u 1:2 w p ls 11 t 'real E',\
'' u 1:3 w p ls 21 t 'real E_{dir}', \
'' u 1:4 w p ls 31 t 'real E_{rec}', \
'esti.spme.ana.error.out' u 1:2 w l ls 11 t 'estimated E', \
'' u 1:3 w l ls 21 t 'estimated E_{dir}', \
'' u 1:4 w l ls 31 t 'estimated E_{rec}'