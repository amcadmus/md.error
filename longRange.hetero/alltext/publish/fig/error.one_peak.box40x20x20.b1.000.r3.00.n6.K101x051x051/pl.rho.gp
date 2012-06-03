set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig.rho.eps'

set style line 11 lc 1 lw 3 pt 7
set style line 21 lc 3 lw 3 pt 7
set style line 31 lc 2 lw 3 pt 7
set style line 41 lc 4 lw 3 pt 7

set format y '%.1f'
set xlabel 'x [ nm ]'
set mxtics 5
#set ylabel '{/Symbol r}(x) [ e/nm^{3} ]'
set ylabel 'Charge densities [ e/nm^{3} ]'
set yrange [-1.2:1.5]
set ytics .5
set mytics 5
set key top right

# pl\
# 'rho.x.avg.out' u 1:2 w l ls 11 title '{/Symbol r}_{+}(x)', \
# '' u 1:($3) w l ls 21 title '{/Symbol r}_{-}(x)',\
# '' u 1:($2+$3) w l ls 31 title '{/Symbol r}_{q}(x)',\
# '' u 1:($2-$3) w l ls 41 title '{/Symbol r}_{q^2}(x)'

pl \
0 ls 0 lw 3 not,\
'rho.x.avg.out' u 1:2 w l ls 11 title '{/Symbol r}_{+}(x)', \
'' u 1:($3) w l ls 21 title '{/Symbol r}_{-}(x)'
