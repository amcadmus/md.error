set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig.rand2.meanf.eps'

set style line 01 lc 0 lw 3 pt 7
set style line 11 lc 1 lw 3 pt 2
set style line 21 lc 2 lw 3 pt 2

set style line 31 lc 3 lw 3 pt 2 lt 1
set style line 32 lc 3 lw 3 pt 5 lt 2
set style line 33 lc 3 lw 3 pt 1
set style line 34 lc 3 lw 3 pt 4 lt 3
set style line 41 lc 4 lw 3 pt 2 lt 1

set format y '%.1e'
set xlabel 'x [ nm ]'
set ylabel 'Mean error force [ kJ/(mol nm) ]'
set xrange [0:14.896480]
set yrange [-12e-2:12e-2]
#set ytics 2e-2
set mxtics 2
set mytics 5

set key bottom right

# set ytics nomirror
# set y2range [-1.2:1.2]
# set y2tics .2
# set my2tics 2
# set format y2 '%.1f'
# set y2label '{/Symbol r}_q(x) [ e/nm^3 ]'
# set rmargin 8.5

set arrow from 5.5, -7e-2 to 4.5, -7e-2
set label "Ewald dir" at 5.7, -7e-2

pl\
'rand2/real.ewald.meanf.out'		u 1:(138.935485*$8) w p ls 21 not 'real Ewald rec',\
'rand2/esti.ewald.meanf.out'		u 1:(138.935485*$8) w l ls 21 not 'estimated Ewald rec',\
'rand2/real.spme.ana.st.meanf.out'	u 1:(138.935485*$5) w p ls 11 not 'real dir',\
'rand2/real.spme.ana.st.meanf.out'	u 1:(138.935485*$8) w p ls 31 not 'real rec ana st', \
'rand2/real.spme.ana.selfCorr.meanf.out'u 1:(138.935485*$8) w p ls 32 not 'real rec ana sc', \
'rand2/real.spme.ik.meanf.out'		u 1:(138.935485*$8) w p ls 34 not 'real rec ik', \
'rand2/real.spme.ik.st.meanf.out'	u 1:(138.935485*$8) w p ls 33 not 'real rec ik st', \
'rand2/esti.spme.ik.meanf.out'		u 1:(138.935485*$8) w l ls 34 not 'esti ik', \
'rand2/esti.spme.st.meanf.out'		u 1:(138.935485*$5) w l ls 11 not 'esti dir',\
'rand2/esti.spme.st.meanf.out'		u 1:(138.935485*$8) w l ls 31 not 'esti rec st', \
'rand2/esti.spme.ana.selfCorr.meanf.out'u 1:(138.935485*$8) w l ls 32 not 'esti rec ana sc'

# 'real.spme.ik.st.meanf.out'		u 1:(138.935485*$8) w p ls 33 not 'real rec ik st', \
# 'real.spme.ana.st.meanf.out'		u 1:(138.935485*$2) w p ls 41 not 'real SPME ana ', \
# 'esti.spme.st.meanf.out'		u 1:(138.935485*$2) w l ls 41 not 'esti SPME ana '


# pl\
# 'real.ewald.meanf.out' u 1:(138.935485*$8) w p ls 21 not 'real Ewald rec',\
# 'esti.ewald.meanf.out' u 1:(138.935485*$8) w l ls 21 not 'esti Ewald rec',\
# 'real.spme.ana.meanf.out' u 1:(138.935485*$5) w p ls 11 not 'real <{/Symbol D}F_x^{dir}>', \
# 'esti.spme.ana.meanf.out' u 1:(138.935485*$5) w l ls 11 not 'estimate <{/Symbol D}F_x^{dir}>', \
# 'real.spme.ana.meanf.out' u 1:(138.935485*$8) w p ls 31 not 'real <{/Symbol D}F_x^{rec}>', \
# 'esti.spme.ana.meanf.out' u 1:(138.935485*$8) w l ls 31 not 'estimate <{/Symbol D}F_x^{rec}>',\
# 'real.spme.ana.meanf.out' u 1:(138.935485*$2) w p ls 41 not 'real <{/Symbol D}F_x>',\
# 'esti.spme.ana.meanf.out' u 1:(138.935485*$2) w l ls 41 not 'estimate <{/Symbol D}F_x>'

# 'rho.x.avg.out' u 1:($2+$3) axes x2y2 w l ls 01 not 'rho',\
