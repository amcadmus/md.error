set term post eps enh color font 16 size 8.5cm,6cm
set out 'fig.ana.ewald.error.eps'

set style line 90 lc 0 lw 3 pt 3 lt 1
set style line 11 lc 1 lw 3 pt 2 lt 1
set style line 21 lc 2 lw 3 pt 2 lt 1
set style line 31 lc 3 lw 3 pt 2 lt 1
set style line 32 lc 3 lw 3 pt 5 lt 2
set style line 33 lc 3 lw 3 pt 1
set style line 34 lc 3 lw 3 pt 7 lt 3
set style line 41 lc 4 lw 3 pt 2 lt 1

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
set mxtics 2
set mytics 2

pl\
'real.ewald.error.out' u 1:(138.935485*$4) w p ls 21 not 'real Ewald rec',\
'esti.ewald.error.out' u 1:(138.935485*$4) w l ls 21 not 'estimated Ewald rec'

set bmargin 0
set tmargin 1
set size 1,.7
set origin 0,.3
set yrange [1e-4:4.5e-2]
set ytics 1e-2
unset xlabel
set format x ""
set mytics 2
set ylabel "RMS error [ kJ/(mol nm) ]"

set arrow from 7.3,1.2e-2 to 6.8,1.5e-2 lw 1
set label "st ik, st ana" at 7.5,1.0e-2
set arrow from 7.3,2.3e-2 to 6.8,2.0e-2 lw 1
set label "orig ik" at 7.5,2.25e-2

pl\
'real.spme.ana.st.error.out'		u 1:(138.935485*$3) w p ls 11 not 'real dir',\
'real.spme.ik.error.out'		u 1:(138.935485*$4) w p ls 34 not 'real rec ik', \
'real.spme.ana.selfCorr.error.out'	u 1:(138.935485*$4) w p ls 32 not 'real rec ana sc', \
'real.spme.ik.st.error.out'		u 1:(138.935485*$4) w p ls 33 not 'real rec ik st', \
'real.spme.ana.st.error.out'		u 1:(138.935485*$4) w p ls 31 not 'real rec ana st', \
'real.spme.ana.st.error.out'		u 1:(138.935485*$2) w p ls 41 not 'real ana ', \
'esti.spme.st.error.out'		u 1:(138.935485*$3) w l ls 11 not 'esti dir',\
'esti.spme.st.error.out'		u 1:(138.935485*$4) w l ls 31 not 'esti rec st', \
'esti.spme.ik.error.out'		u 1:(138.935485*$4) w l ls 34 not 'esti ik', \
'esti.spme.ana.selfCorr.error.out'	u 1:(138.935485*$4) w l ls 32 not 'esti rec ana sc', \
'esti.spme.st.error.out'		u 1:(138.935485*$2) w l ls 41 not 'esti st'

unset multiplot


#'real.spme.ik.error.out'		u 1:(138.935485*$4) w p ls 34 t 'real SPME ik rec', \
