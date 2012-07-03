set term post eps enh color  font 16 size 8.5cm,6cm
set out 'fig.rho.eps'

set style line 11 lc 1 lw 3 pt 7 lt 1
set style line 12 lc 1 lw 3 pt 7 lt 2
set style line 21 lc 3 lw 3 pt 7 lt 1

set format y '%.1f'
set xlabel 'x [ nm ]'
set mxtics 2
set xrange [0:14.896480]
set ylabel 'Charge number density [ nm^{-3} ]'
set yrange [-.1:90]
#set ytics .5
set mytics 5
#set key at 13,70 cent

pl \
'rand1/rho.x.avg.out' u 1:($2/0.834) w l ls 21 title '-',\
'rand1/rho.x.avg.out' u 1:(-$3/0.417) w l ls 11 title 'Example 1  +',\
'rand2/rho.x.avg.out' u 1:(-$3/0.417) w l ls 12 title 'Example 2  +'

#0 ls 0 lw 3 not,\
