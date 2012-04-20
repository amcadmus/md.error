set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig.ana.ewald.error.eps'

set style line 11 lc 1 lw 3 pt 2
set style line 21 lc 2 lw 3 pt 2
set style line 31 lc 3 lw 3 pt 2
set style line 32 lc 3 lw 3 pt 5 lt 0
set style line 33 lc 3 lw 3 pt 1
set style line 34 lc 3 lw 3 pt 7
set style line 41 lc 4 lw 3 pt 2

set multiplot
set size 1,.3
set origin 0,0.0

set bmargin 3
set tmargin 0
set lmargin 12
set rmargin 2

set format y '%.1e'
set xlabel 'x [ nm ]'
set yrange [0:6.0e-10]
unset ylabel
set ytics 2e-10
set mxtics 5
set mytics 2

pl\
'real.ewald.error.out' u 1:(138.935485*$4) w p ls 21 not 'real Ewald rec',\
'esti.ewald.error.out' u 1:(138.935485*$4) w l ls 21 not 'estimated Ewald rec'

set bmargin 0
set tmargin 1
set size 1,.7
set origin 0,.3
set yrange [4e-4:14e-2]
set ytics 2e-2
unset xlabel
set format x ""
set mytics 2
set ylabel "RMS error [ kJ/(mol nm) ]"

pl\
'real.spme.ana.st.error.out'		u 1:(138.935485*$3) w p ls 11 title 'real dir',\
'real.spme.ana.st.error.out'		u 1:(138.935485*$4) w p ls 31 title 'real rec ana st', \
'real.spme.ik.st.error.out'		u 1:(138.935485*$4) w p ls 33 title 'real rec ik st', \
'real.spme.ana.selfCorr.error.out'	u 1:(138.935485*$4) w p ls 32 title 'real rec ana sc', \
'real.spme.ana.st.error.out'		u 1:(138.935485*$2) w p ls 41 title 'real ana ', \
'esti.spme.st.error.out'		u 1:(138.935485*$3) w l ls 11 title 'esti dir',\
'esti.spme.st.error.out'		u 1:(138.935485*$4) w l ls 31 title 'esti rec st', \
'esti.spme.ana.selfCorr.error.out'	u 1:(138.935485*$4) w l ls 32 title 'esti rec ana sc', \
'esti.spme.st.error.out'		u 1:(138.935485*$2) w l ls 41 title 'esti st'

unset multiplot


#'real.spme.ik.error.out'		u 1:(138.935485*$4) w p ls 34 t 'real SPME ik rec', \
