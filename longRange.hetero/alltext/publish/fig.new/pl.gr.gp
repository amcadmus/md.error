set term post eps enh color  font 16 size 8.5cm,6cm
set out 'fig.gr.eps'

set style line 11 lc 1 lw 3 pt 7 lt 1
set style line 12 lc 1 lw 3 pt 7 lt 2
set style line 21 lc 0 lw 3 pt 7 lt 1

set format y '%.1e'
set xlabel 'R_g [ nm ]'
set mxtics 2
#set xrange [0:14.896480]
set ylabel 'RMS error [ kJ/(mol nm) ]'
# #set ytics .5
# set mytics 5
# #set key at 13,70 cent

pl \
'gr.out' u 1:(138.935485*$2) w l ls 21 notitle '-'

#0 ls 0 lw 3 not,\
