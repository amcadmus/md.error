set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig.water.ik.error.eps'

set style line 11 lc 1 lw 3 pt 1 lt 1
set style line 21 lc 2 lw 3 pt 1 lt 2
set style line 31 lc 3 lw 3 pt 1 lt 3
set style line 41 lc 3 lw 3 pt 1 lt 0
set style line 51 lc 4 lw 3 pt 1 lt 0
set style line 12 lc 1 lw 3 pt 2
set style line 22 lc 2 lw 3 pt 2
set style line 32 lc 3 lw 3 pt 2
set style line 13 lc 1 lw 3 pt 3
set style line 23 lc 2 lw 3 pt 3
set style line 33 lc 3 lw 3 pt 3
set style line 14 lc 1 lw 3 pt 4
set style line 24 lc 2 lw 3 pt 4
set style line 34 lc 3 lw 3 pt 4

set format y '%.1e'
set xlabel 'x'
set ylabel 'RMS error'
#set logscale y

pl\
'real.spme.ik.error.out' 	u 1:4 w p ls 31 t 'real E',\
'real.spme.ik.st.error.out'	u 1:4 w p ls 32 t 'st real E',\
'real.spme.ik.noHO.error.out'	u 1:4 w p ls 33 t 'noHO real E',\
'real.spme.ik.noHO.st.error.out'u 1:2 w p ls 14 t 'noHO st real E',\
'real.spme.ik.noHO.st.error.out'u 1:3 w p ls 24 t 'noHO st real E',\
'real.spme.ik.noHO.st.error.out'u 1:4 w p ls 34 t 'noHO st real E',\
'esti.spme.ik.error.out' u 1:4 w l ls 51 t 'estimated E_{rec}',\
'esti.spme.st.error.out' u 1:4 w l ls 41 t 'estimated st E_{rec}',\
'esti.spme.st.h2o.error.out' u 1:2 w l ls 11 t 'estimated st h2o E', \
'esti.spme.st.h2o.error.out' u 1:3 w l ls 21 t 'estimated st h2o E_{dir}', \
'esti.spme.st.h2o.error.out' u 1:4 w l ls 31 t 'estimated st h2o E_{rec}'

