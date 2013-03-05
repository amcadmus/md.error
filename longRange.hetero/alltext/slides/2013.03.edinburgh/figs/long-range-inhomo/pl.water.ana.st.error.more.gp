set term post eps enh color font 16 size 8.5cm,6cm
set out 'fig.water.ana.st.error.more.eps'

set style line 90 lc 0 lw 3 pt 3 lt 1
set style line 11 lc 1 lw 3 pt 2 lt 1
set style line 21 lc 2 lw 3 pt 1 lt 1
set style line 31 lc 3 lw 3 pt 2 lt 1
set style line 32 lc 3 lw 3 pt 3 lt 2
set style line 33 lc 3 lw 3 pt 1
set style line 34 lc 3 lw 3 pt 4 lt 3
set style line 35 lc 3 lw 3 pt 7 lt 3
set style line 41 lc 4 lw 3 pt 2 lt 1

set format y '%.1e'
set xlabel 'x [ nm ]'
set ylabel "RMS error [ kJ/(mol nm) ]"
set xrange [0:14.896480]
set yrange [0:1.8e-2]
# set ytics 5e-3
set mxtics 2
set mytics 2

set arrow from 6,1.55e-2 to 6,0.45e-2 lw 1
set arrow from 7,0.43e-2 to 7,0.18e-2 lw 3 lt 3 lc 0
# set label "1st nbr. approx." at 6.2, 0.6e-2
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
'water/real.spme.ana.noHO.st.error.out'	u 1:(138.935485*$4) w p ls 31 not 'real rec ana st', \
'water/esti.spme.st.h2o.error.out'	u 1:(138.935485*$4) w l ls 31 not 'esti rec st', \
'water/esti.spme.st.h2o.gr1.8.error.out'	u 1:(138.935485*$4) w l ls 34 not 'esti rec st gr 1.8', \
'water/esti.spme.st.error.out'		u 1:(138.935485*$4) w l ls 32 not 'esti rec ana sc',\
'water/real.spme.ik.noHO.st.error.out'	u 1:(138.935485*$4) w p ls 21 not 'real rec ik st'


unset multiplot


