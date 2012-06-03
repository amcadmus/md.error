set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig.ik.meanf.eps'

set style line 11 lc 1 lw 3 pt 7
set style line 21 lc 2 lw 3 pt 7
set style line 31 lc 3 lw 3 pt 7

set format y '%.1e'
set xlabel 'x'
set ylabel 'Mean error force'
set yrange [-.5e-3:.5e-3]
set ytics 1e-4
set mxtics 5
set mytics 5

set key bottom right

pl\
'real.spme.ik.meanf.out' u 1:2 w p ls 11 t 'real <{/Symbol D}F_x>',\
'' u 1:5 w p ls 21 t 'real <{/Symbol D}F_x^{dir}>', \
'' u 1:8 w p ls 31 t 'real <{/Symbol D}F_x^{rec}>', \
'esti.spme.ik.meanf.out' u 1:2 w l ls 11 t 'estimate <{/Symbol D}F_x>',\
'' u 1:5 w l ls 21 t 'estimate <{/Symbol D}F_x^{dir}>', \
'' u 1:8 w l ls 31 t 'estimate <{/Symbol D}F_x^{rec}>'
