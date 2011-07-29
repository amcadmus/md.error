#set border 4095 front linetype -1 linewidth 1.000
set term post eps enh color solid font 16
set out 'error-rcut.eps'

# unset xtics
# unset ytics
# set tmargin 0n
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
set pm3d implicit at b
#set palette defined (2.5 "blue", 3.75 "white", 5 "red")
set palette rgbformulae 22,13,-31
#set size ratio 0.88
#set size ratio 0.90

set size 2,1.36
set origin 0, 0
set multiplot 
set size 0.5,0.72

set cbrange [0:0.035]
set format cb "%.3f"

set origin 0, 0.68
set title 't^*=20' font "Helvetica,24"
spl'dynamic.uniform.rc05.0/error.t0020.00.out' not

set origin 0.5, 0.68
set title 't^*=32' font "Helvetica,24"
spl'dynamic.uniform.rc05.0/error.t0032.00.out' not

set origin 1.0, 0.68
set title 't^*=48'  font "Helvetica,24"
spl'dynamic.uniform.rc05.0/error.t0048.00.out' not

set origin 1.5, 0.68
set title 't^*=80' font "Helvetica,24"
spl'dynamic.uniform.rc05.0/error.t0080.00.out' not

set cbrange [2.5:5]
set format cb "%.1f"

set origin 0.0, 0.0
set title 't^*=20' font "Helvetica,24"
spl'dynamic.adapt-new/rcut.t0020.00.out' not

set origin 0.5, 0.0
set title 't^*=32' font "Helvetica,24"
spl'dynamic.adapt-new/rcut.t0032.00.out' not

set origin 1.0, 0.0
set title 't^*=48' font "Helvetica,24"
spl'dynamic.adapt-new/rcut.t0048.00.out' not

set origin 1.5, 0.0
set title 't*=80' font "Helvetica,24"
spl'dynamic.adapt-new/rcut.t0080.00.out' not


unset multiplot
