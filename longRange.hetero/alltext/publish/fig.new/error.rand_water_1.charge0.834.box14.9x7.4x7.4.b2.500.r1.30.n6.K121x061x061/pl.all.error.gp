set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig.all.error.eps'

set style line 11 lc 1 lw 3 pt 7
set style line 21 lc 2 lw 3 pt 7
set style line 31 lc 3 lw 3 pt 7
set style line 41 lc 4 lw 3 pt 7

set multiplot
set size 1,.4
set origin 0,0.0

set format y '%.0e'
set xlabel 'x [ nm ]'
set yrange [0:4.0e-8]
set ytics 2e-8
set mxtics 5
set mytics 4

pl\
'real.ewald.error.out' u 1:4 w p ls 21 not 'real Ewald rec',\
'esti.ewald.error.out' u 1:4 w l ls 21 not 'estimated Ewald rec'

set size 1,.7
set origin 0,.295
set yrange [4e-6:5e-4]
set ytics 1e-4
unset xlabel
set format x ""
set mytics 5

pl\
'real.ewald.error.out' u 1:3 w p ls 11 not 'real dir',\
'esti.ewald.error.out' u 1:3 w l ls 11 not 'estimated dir',\
'real.spme.ik.error.out' u 1:4 w p ls 31 not 'real SPME ik rec', \
'esti.spme.ik.error.out' u 1:4 w l ls 31 not 'estimated SPME ik rec', \
'real.spme.ana.error.out' u 1:4 w p ls 41 not 'real SPME ana rec', \
'esti.spme.ana.error.out' u 1:4 w l ls 41 not 'estimated SPME ana rec'
