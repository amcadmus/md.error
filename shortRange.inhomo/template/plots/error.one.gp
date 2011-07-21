#set border 4095 front linetype -1 linewidth 1.000
set view map
set samples 50, 50
set isosamples 50, 50
unset surface
set style data pm3d
set style function pm3d
set ticslevel 0
set xlabel "x" 
set xrange [ 0 : 150 ] noreverse nowriteback
set ylabel "y" 
set yrange [ 0:140 ] noreverse nowriteback
set pm3d implicit at b
set size ratio 0.88
set zrange [0:006]

set term gif
set out ''
spl '' not
