set term post eps enh solid  font 16 size 8.5cm,6cm
set out 'natom-rcut-c.eps'

set style line 1 lc 1 lw 2 pt 7
set style line 2 lc 3 lw 2 pt 7
set style line 3 lc 4 lw 2 pt 7

set mxtics 2
set mytics 2
set xrange[3:12]
#set yrange[32:40]
set xlabel 'r_c^*'
set ylabel 'c_P^*'

pl \
'rcut-N00004000-P0.1535.out' u 1:2:($3*2) w e ls 1 t 'N=04,000',\
'' w l ls 1 not, \
'rcut-N00016000-P0.1535.out' u 1:2:($3*2) w e ls 2 t 'N=16,000',\
'' w l ls 2 not,\
'rcut-N00024000-P0.1535.out' u 1:2:($3*2) w e ls 3 t 'N=24,000',\
'' w l ls 3 not


