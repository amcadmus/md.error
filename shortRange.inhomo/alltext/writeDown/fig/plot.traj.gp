set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'trajs.eps'

set style line 1 lt 1 lc 1 lw 3 pt 5
set style line 2 lt 1 lc 2 lw 3 pt 5
set style line 3 lt 1 lc 3 lw 3 pt 5
set style line 4 lt 1 lc 4 lw 3 pt 5
set style line 5 lt 1 lc 5 lw 3 pt 5


set xrange [40:80]
set yrange [50:77]
set xlabel 'x'
set ylabel 'y'
set mxtics 5
set mytics 5
set format y "%.0f"
set key bottom left

pl\
'collision-x0.60-u2.20-rc2.5/traj.a.out' u 1:3 w lp ls 1 t 'URC 2.5',\
'collision-x0.60-u2.20-rc5.0/traj.a.out' u 1:3 w lp ls 2 t 'URC 5.0',\
'collision-x0.60-u2.20-rc7.5/traj.a.out' u 1:3 w lp ls 3 t 'URC 7.5',\
'adaptColl-x0.60-u2.20-rc5.0-61/traj.a.out' u 1:3 w lp ls 4 t 'ARC 5.0',\
'fcorrColl-x0.60-u2.20-rc2.5-14/traj.a.out' u 1:3 w lp ls 5 t 'LFC 2.5'