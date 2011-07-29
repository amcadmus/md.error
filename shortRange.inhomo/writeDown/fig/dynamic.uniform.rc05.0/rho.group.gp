#set border 4095 front linetype -1 linewidth 1.000
set term post eps enh color solid font 16
set out 'rho.group.eps'

# unset xtics
# unset ytics
# set tmargin 0
# set bmargin 0
# set lmargin 2
# set rmargin 2
# set multiplot layout 2,3
#set multiplot layout 1,2 title "Auto-layout of stacked plots\n"
unset title

set view map
set samples 50, 50
set isosamples 50, 50
unset surface
set style data pm3d
set style function pm3d
set ticslevel 0
set xrange [ 20 : 130 ] noreverse nowriteback
set yrange [ 20:120 ] noreverse nowriteback
set cbrange [0:1.2]
set pm3d implicit at b
#set palette defined (2.5 "blue", 3.75 "white", 5 "red")
set palette rgbformulae 22,13,-31
#set size ratio 0.88
#set size ratio 0.90

set size 2,1.36
set origin 0, 0
set multiplot 
set size 0.5,0.72

set origin 0, 0.68
set title 't^*=20' font "Helvetica,24"
spl'rho.t0020.00.out' not

set origin 0.5, 0.68
set title 't^*=32' font "Helvetica,24"
spl'rho.t0032.00.out' not
# set xtics nomirror
# set tics scale 0

set origin 1.0, 0.68
set title 't^*=48'  font "Helvetica,24"
spl'rho.t0048.00.out' not

set origin 1.5, 0.68
set title 't^*=64' font "Helvetica,24"
spl'rho.t0064.00.out' not

set origin 0.0, 0.0
set title 't^*=80' font "Helvetica,24"
spl'rho.t0080.00.out' not

set origin 0.5, 0.0
set title 't^*=100' font "Helvetica,24"
spl'rho.t0100.00.out' not

set origin 1.0, 0.0
set title 't^*=120' font "Helvetica,24"
spl'rho.t0120.00.out' not

set origin 1.5, 0.0
set title 't*=140' font "Helvetica,24"
spl'rho.t0140.00.out' not

unset multiplot
