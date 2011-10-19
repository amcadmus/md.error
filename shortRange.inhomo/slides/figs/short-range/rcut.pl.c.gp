set term post eps enh solid  font 16 size 6cm,4.5cm
set out 'natom-rcut-c.eps'

set style line 1 lc 1 lw 2 pt 7
set style line 2 lc 3 lw 2 pt 7

set xtics 1
set mxtics 1
set mytics 5
set xrange[2:12]
set yrange[45:70]
set xlabel 'r_c^*'
set ylabel 'c_P^*'

pl \
'rcut-N00016000-P0.1535.out' u 1:2:($3*2) w e ls 1 not,\
'' w l ls 1 not


