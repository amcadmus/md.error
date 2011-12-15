#set border 4095 front linetype -1 linewidth 1.000
set view map
set samples 50, 50
set isosamples 50, 50
unset surface
set style data pm3d
set style function pm3d
set ticslevel 0
set xlabel "x" 
set xrange [ 18 : 110 ] noreverse nowriteback
set ylabel "y" 
set yrange [ 18 : 110 ] noreverse nowriteback
set pm3d implicit at b
#set size ratio 0.88
set size ratio 1.0

set term gif
set out ''
spl '' not
