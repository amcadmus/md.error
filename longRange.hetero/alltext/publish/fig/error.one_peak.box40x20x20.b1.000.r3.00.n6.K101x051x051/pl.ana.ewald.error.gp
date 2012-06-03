set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig.ana.ewald.error.eps'

set style line 11 lc 1 lw 3 pt 7
set style line 21 lc 2 lw 3 pt 7
set style line 31 lc 3 lw 3 pt 7
set style line 41 lc 4 lw 3 pt 7

set multiplot
set size 1,.3
set origin 0,0.0

set bmargin 3
set tmargin 0
set lmargin 10
set rmargin 2

set format y '%.0e'
set xlabel 'x [ nm ]'
set yrange [0:6.0e-6]
unset ylabel
set ytics 2e-6
set mxtics 5
set mytics 2

pl\
'real.ewald.error.out' u 1:(138.935485*$4) w p ls 21 not 'real Ewald rec',\
'esti.ewald.error.out' u 1:(138.935485*$4) w l ls 21 not 'estimated Ewald rec'

set bmargin 0
set tmargin 1
set size 1,.7
set origin 0,.3
set yrange [4e-4:6e-2]
set ytics 1e-2
unset xlabel
set format x ""
set mytics 5
set ylabel "RMS error [ kJ/(mol nm) ]"

pl\
'real.ewald.error.out' u 1:(138.935485*$3) w p ls 11 not 'real dir',\
'esti.ewald.error.out' u 1:(138.935485*$3) w l ls 11 not 'estimated dir',\
'real.spme.ana.error.out' u 1:(138.935485*$4) w p ls 31 not 'real SPME ana rec', \
'esti.spme.ana.error.out' u 1:(138.935485*$4) w l ls 31 not 'estimated SPME ik rec', \
'real.spme.ana.error.out' u 1:(138.935485*$2) w p ls 41 not 'real SPME ana ', \
'esti.spme.ana.error.out' u 1:(138.935485*$2) w l ls 41 not 'estimated SPME ana '

unset multiplot
