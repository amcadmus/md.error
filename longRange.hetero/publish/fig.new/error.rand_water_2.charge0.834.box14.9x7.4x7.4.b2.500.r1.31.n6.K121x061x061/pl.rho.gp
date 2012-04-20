set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig.rho.eps'

set style line 11 lc 1 lw 3 pt 7
set style line 21 lc 3 lw 3 pt 7

set format y '%.1f'
set xlabel 'x [ nm ]'
set mxtics 2
set ylabel 'Charge number density [ nm^{-3} ]'
set xrange [0:14.896480]
set yrange [-.1:80]
#set ytics .5
set mytics 5
set key bottom right
set key at 13.2,30 cent

pl \
'rho.x.avg.out' u 1:($2/0.834) w l ls 21 title '-',\
'rho.x.avg.out' u 1:(-$3/0.417) w l ls 11 title '+'

#0 ls 0 lw 3 not,\
