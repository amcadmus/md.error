set term post eps enh color font 12 size 8.5cm,6.5cm
set out 'fig.order.st.eps'

set style line 11 lc 1 lw 3 pt 2 lt 1
set style line 12 lc 1 lw 3 pt 2 lt 2
set style line 21 lc 2 lw 3 pt 1 lt 1
set style line 22 lc 2 lw 3 pt 2 lt 2
set style line 25 lc 2 lw 3 pt 4 lt 1
set style line 31 lc 3 lw 3 pt 2 lt 1
set style line 32 lc 3 lw 3 pt 2 lt 2
set style line 33 lc 3 lw 3 pt 1 lt 1
set style line 35 lc 3 lw 3 pt 4 lt 1

set format x '%.1f'
set format y '10^{%L}'
set logscale y

set mxtics 5
set ytics 10
set mytics 10
set xrange [.9:5.1]
set yrange [1e-7:3e2]
set xlabel '{/Symbol b} [ nm^{-1} ]'
set ylabel "RMS error [ kJ/(mol nm) ]"

# 'r1.300.n8.K041x011x011.st.noHO.error' u 1:(138.935485*$7) ls 31 w p not,\
# '' u 1:(138.935485*$4)  ls 21 w p not,\
# '' u 1:(138.935485*$10) ls 32 w l not,\
# '' u 1:(138.935485*$13) ls 31 w l not,\

# 'r1.300.n8.K081x021x021.st.noHO.error' u 1:(138.935485*$7) ls 31 w p not,\
# '' u 1:(138.935485*$4)  ls 21 w p not,\
# '' u 1:(138.935485*$10) ls 32 w l not,\
# '' u 1:(138.935485*$13) ls 31 w l not,\
# 'r1.300.n8.K121x031x031.st.noHO.error' u 1:(138.935485*$7) ls 31 w p not,\
# '' u 1:(138.935485*$4)  ls 21 w p not,\
# '' u 1:(138.935485*$10) ls 32 w l not,\
# '' u 1:(138.935485*$13) ls 31 w l not,\

set label "0.9 nm" at 4.5,5e-6
set label "1.1 nm" at 4.3,3e-7
set label "1.3 nm" at 3.65,3e-7
set label "r_c = 1.7 nm" at 2.8,3e-7
set label "n = 4" at 1.0,6e-3
set label "6" at 1.2,3e-5
set label "8" at 1.2,3e-7


pl\
'r0.900.n4.K060x015x015.ik.noHO.error' u 1:(138.935485*$3) ls 12 w p not,\
'' u 1:(138.935485*$6)  ls 12 w l not,\
'r1.100.n4.K060x015x015.ik.noHO.error' u 1:(138.935485*$3) ls 12 w p not,\
'' u 1:(138.935485*$6)  ls 12 w l not,\
'r1.300.n4.K060x015x015.ik.noHO.error' u 1:(138.935485*$3) ls 12 w p not,\
'' u 1:(138.935485*$6)  ls 12 w l not,\
'r1.700.n4.K060x015x015.ik.noHO.error' u 1:(138.935485*$3) ls 12 w p not,\
'' u 1:(138.935485*$6)  ls 12 w l not,\
'rand.r1.300.n4.K121x031x031.st.noHO.error' u 1:(138.935485*$7) ls 35 w p not,\
'rand.r1.300.n6.K121x031x031.st.noHO.error' u 1:(138.935485*$7) ls 35 w p not,\
'rand.r1.300.n8.K121x031x031.st.noHO.error' u 1:(138.935485*$7) ls 35 w p not,\
'r1.300.n4.K120x030x030.st.noHO.error' u 1:(138.935485*$7) ls 31 w p not,\
'' u 1:(138.935485*$4)  ls 21 w p not,\
'' u 1:(138.935485*$10) ls 32 w l not,\
'' u 1:(138.935485*$13) ls 31 w l not,\
'r1.300.n6.K120x030x030.st.noHO.error' u 1:(138.935485*$7) ls 31 w p not,\
'' u 1:(138.935485*$4)  ls 21 w p not,\
'' u 1:(138.935485*$10) ls 32 w l not,\
'' u 1:(138.935485*$13) ls 31 w l not,\
'r1.300.n8.K120x030x030.st.noHO.error' u 1:(138.935485*$7) ls 31 w p not,\
'' u 1:(138.935485*$4)  ls 21 w p not,\
'' u 1:(138.935485*$10) ls 32 w l not,\
'' u 1:(138.935485*$13) ls 31 w l not

# '' u 1:(138.935485*$4)  ls 35 w p not,\

# 'rand.r1.300.n6.K120x030x030.st.noHO.error' u 1:(138.935485*$7) ls 25 w p not,\
# '' u 1:(138.935485*$4)  ls 25 w p not,\


# 'r1.300.n8.K081x021x021.ana.noHO.error' u 1:(138.935485*$4) ls 31 w p not,\
# '' u 1:(138.935485*$10) ls 31 w l not,\
# '' u 1:(138.935485*$7)  ls 32 w l not
