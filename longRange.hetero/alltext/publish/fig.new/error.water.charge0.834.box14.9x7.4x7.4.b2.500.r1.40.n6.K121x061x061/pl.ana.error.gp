set term post eps enh color font 16 size 8.5cm,6cm
set out 'fig.ana.error.eps'

set style line 90 lc 0 lw 3 pt 3 lt 1
set style line 11 lc 1 lw 3 pt 2 lt 1
set style line 21 lc 2 lw 3 pt 2 lt 1
set style line 31 lc 3 lw 3 pt 2 lt 1
set style line 32 lc 3 lw 3 pt 5 lt 2
set style line 33 lc 3 lw 3 pt 1
set style line 34 lc 3 lw 3 pt 7 lt 3
set style line 41 lc 4 lw 3 pt 2 lt 1

set format y '%.1e'
set xlabel 'x [ nm ]'
set ylabel "RMS error [ kJ/(mol nm) ]"
set xrange [0:14.896480]
set yrange [0:2e-2]
set ytics 5e-3
set mytics 5

# set arrow from 7.3,0.7e-2 to 6.8,1.0e-2 lw 1
# set label "st ana" at 7.5,0.7e-2
# # set label "st ana" at 7.5,0.8e-2
# set arrow from 7.3,4.2e-2 to 6.8,3.4e-2 lw 1
# set label "orig ana" at 7.5,4.5e-2
# set arrow from 5.5, 7e-2 to 4.5, 7e-2
# set label "Ewald dir" at 5.7, 7e-2
# # set arrow from 10.6,3.9e-2 to 10.3,1.8e-2 lw .5
# # set label "orig ik" at 10.5,4.3e-2
# # set label "Ewald dir" at 7.5,3.1e-2
# # set arrow from 7.3,4.1e-2 to 6.8,3.8e-2 lw 1
# # set label "orig ana" at 7.5,4.1e-2

pl\
'real.spme.ana.noHO.st.error.out'	u 1:(138.935485*$3) w p ls 11 not 'real dir',\
'real.spme.ana.noHO.st.error.out'	u 1:(138.935485*$4) w p ls 31 not 'real rec ana st', \
'esti.spme.st.h2o.error.out'		u 1:(138.935485*$3) w l ls 11 not 'esti dir',\
'esti.spme.st.h2o.error.out'		u 1:(138.935485*$4) w l ls 31 not 'esti rec st', \
'esti.spme.st.error.out'		u 1:(138.935485*$4) w l ls 32 not 'esti rec ana sc',\
'tmp.out'	u 1:(138.935485*$4) w p lc 0 lw 3 not 'esti rec ana sc'

# 'esti.spme.ana.selfCorr.error.out'	u 1:(138.935485*$4) w l ls 32 not 'esti rec ana sc',\
# 'real.spme.ana.noHO.error.out'		u 1:(138.935485*$4) w p ls 32 not 'real rec ana sc', \
# 'real.spme.ik.st.error.out'		u 1:(138.935485*$4) w p ls 33 not 'real rec ik st', \
# 'real.spme.ik.error.out'		u 1:(138.935485*$4) w p ls 34 not 'real rec ik', \
# 'esti.spme.ik.error.out'		u 1:(138.935485*$4) w l ls 34 not 'esti ik', \

# 'real.spme.ana.st.error.out'		u 1:(138.935485*$2) w p ls 41 not 'real ana ', \
# 'esti.spme.st.error.out'		u 1:(138.935485*$2) w l ls 41 not 'esti st'

unset multiplot


#'real.spme.ik.error.out'		u 1:(138.935485*$4) w p ls 34 t 'real SPME ik rec', \
