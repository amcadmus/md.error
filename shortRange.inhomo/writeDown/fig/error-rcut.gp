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
set xrange [ 18 : 110 ] noreverse nowriteback
set xtics 20
set mxtics 2
set yrange [ 18 : 110 ] noreverse nowriteback
set ytics 20
set mytics 2
set pm3d implicit at b
#set palette defined (2.5 "blue", 3.75 "white", 5 "red")
set palette rgbformulae 22,13,-31
#set size ratio 0.88
#set size ratio 0.90

set size 2.05,1.40
set origin 0, 0
set multiplot 
set size 0.5,0.78

set cbrange [0:0.035]
set format cb "%.3f"

# old: 12 32 48 80
set origin 0, 0.68
set title 't^*=12' font "Helvetica,24"
spl'collision-x0.60-u2.20-rc5.0/plots/error.t0012.00.out' not

set origin 0.5, 0.68
set title 't^*=24' font "Helvetica,24"
spl'collision-x0.60-u2.20-rc5.0/plots/error.t0024.00.out' not

set origin 1.0, 0.68
set title 't^*=36'  font "Helvetica,24"
spl'collision-x0.60-u2.20-rc5.0/plots/error.t0036.00.out' not

set origin 1.5, 0.68
set title 't^*=56' font "Helvetica,24"
spl'collision-x0.60-u2.20-rc5.0/plots/error.t0056.00.out' not

set cbrange [2.5:5]
set format cb "%.1f"

set origin 0.0, 0.0
set title 't^*=12' font "Helvetica,24"
spl'adaptColl-x0.60-u2.20-rc5.0-61/plots/rcut.t0012.00.out' not

set origin 0.5, 0.0
set title 't^*=24' font "Helvetica,24"
spl'adaptColl-x0.60-u2.20-rc5.0-61/plots/rcut.t0024.00.out' not

set origin 1.0, 0.0
set title 't^*=36' font "Helvetica,24"
spl'adaptColl-x0.60-u2.20-rc5.0-61/plots/rcut.t0036.00.out' not

set origin 1.5, 0.0
set title 't*=56' font "Helvetica,24"
spl'adaptColl-x0.60-u2.20-rc5.0-61/plots/rcut.t0056.00.out' not


unset multiplot
