set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig.ana.ewald.meanf.eps'

set style line 01 lc 0 lw 3 pt 7
set style line 11 lc 1 lw 3 pt 7
set style line 21 lc 2 lw 3 pt 7
set style line 31 lc 3 lw 3 pt 7
set style line 41 lc 4 lw 3 pt 7

set format y '%.0e'
set xlabel 'x [ nm ]'
set ylabel 'Mean error force [ kJ/(mol nm) ]'
set yrange [-6e-2:6e-2]
set ytics 2e-2
set mxtics 5
set mytics 5

set key bottom right

# set ytics nomirror
# set y2range [-1.2:1.2]
# set y2tics .2
# set my2tics 2
# set format y2 '%.1f'
# set y2label '{/Symbol r}_q(x) [ e/nm^3 ]'
# set rmargin 8.5

pl\
'real.ewald.meanf.out' u 1:(138.935485*$8) w p ls 21 not 'real Ewald rec',\
'esti.ewald.meanf.out' u 1:(138.935485*$8) w l ls 21 not 'estimated Ewald rec',\
'real.spme.ana.meanf.out' u 1:(138.935485*$5) w p ls 11 not 'real <{/Symbol D}F_x^{dir}>', \
'esti.spme.ana.meanf.out' u 1:(138.935485*$5) w l ls 11 not 'estimate <{/Symbol D}F_x^{dir}>', \
'real.spme.ana.meanf.out' u 1:(138.935485*$8) w p ls 31 not 'real <{/Symbol D}F_x^{rec}>', \
'esti.spme.ana.meanf.out' u 1:(138.935485*$8) w l ls 31 not 'estimate <{/Symbol D}F_x^{rec}>',\
'real.spme.ana.meanf.out' u 1:(138.935485*$2) w p ls 41 not 'real <{/Symbol D}F_x>',\
'esti.spme.ana.meanf.out' u 1:(138.935485*$2) w l ls 41 not 'estimate <{/Symbol D}F_x>'

# 'rho.x.avg.out' u 1:($2+$3) axes x2y2 w l ls 01 not 'rho',\
