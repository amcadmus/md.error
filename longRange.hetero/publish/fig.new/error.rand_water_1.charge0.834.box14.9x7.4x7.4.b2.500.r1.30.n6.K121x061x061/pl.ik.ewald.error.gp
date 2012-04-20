set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig.ik.ewald.error.eps'

set style line 11 lc 1 lw 3 pt 7
set style line 21 lc 2 lw 3 pt 7
set style line 31 lc 3 lw 3 pt 7
set style line 32 lc 3 lw 3 pt 2 lt 0
set style line 41 lc 4 lw 3 pt 7

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
set mytics 4
set ylabel "RMS error [ kJ/(mol nm) ]"

pl\
'real.spme.ik.st.error.out' u 1:(138.935485*$3) w p ls 11 not 'real dir',\
'real.spme.ik.st.error.out' u 1:(138.935485*$4) w p ls 31 not 'real SPME ik rec', \
'real.spme.ik.st.error.out' u 1:(138.935485*$2) w p ls 41 not 'real SPME ik ', \
'esti.spme.st.error.out' u 1:(138.935485*$3) w l ls 11 not 'estimated dir',\
'esti.spme.st.error.out' u 1:(138.935485*$4) w l ls 31 not 'estimated SPME ik rec', \
'esti.spme.st.error.out' u 1:(138.935485*$2) w l ls 41 not 'estimated SPME ik '

unset multiplot
