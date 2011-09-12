set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'dists.eps'

set style line 1 lt 1 lc 1 lw 3 pt 5
set style line 2 lt 1 lc 2 lw 3 pt 5
set style line 3 lt 1 lc 3 lw 3 pt 5
set style line 4 lt 1 lc 4 lw 3 pt 5
set style line 5 lt 1 lc 5 lw 3 pt 5


set xrange [0:130]
set yrange [20:30]
set mxtics 2
set mytics 2
set format y "%.0f"
# set key bottom left
set xlabel 't*'
set ylabel 'Distance of COMs'
set key at 125,27

pl\
'collision-x0.60-u2.20-rc2.5/dist.cm.out' u 1:2 w l ls 1 t 'URC 2.5',\
'collision-x0.60-u2.20-rc5.0/dist.cm.out' u 1:2 w l ls 2 t 'URC 5.0',\
'collision-x0.60-u2.20-rc7.5/dist.cm.out' u 1:2 w l ls 3 t 'URC 7.5',\
'adaptColl-x0.60-u2.20-rc5.0-61/dist.cm.out' u 1:2 w l ls 4 t 'ARC 5.0',\
'fcorrColl-x0.60-u2.20-rc2.5-14/dist.cm.out' u 1:2 w l ls 5 t 'FR 2.5'
