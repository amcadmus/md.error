#set border 4095 front linetype -1 linewidth 1.000
set view map
set samples 50, 50
set isosamples 50, 50
unset surface
set style data pm3d
set style function pm3d
set ticslevel 0
set xlabel "x" 
set xrange [ 20 : 130 ] noreverse nowriteback
set ylabel "y" 
set yrange [ 20:120 ] noreverse nowriteback
set cbrange [0:1.2]
set pm3d implicit at b
#set size ratio 0.88
set size ratio 0.90

set term gif
set out ''
spl '' not
