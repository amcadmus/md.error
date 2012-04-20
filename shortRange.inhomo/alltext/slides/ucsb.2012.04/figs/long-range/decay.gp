set term post eps color enh solid
set out 'decay.eps'
set size 0.5,0.5
pl[0:10][0:3] 1/x w l lw 3 title '1/r', 1/(x*x*x*x*x*x) w l lw 3 title '1/r^6'

