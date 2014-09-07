set term post enh eps color solid size 25.5cm,18cm font "Helvetica, 40"
set out 'rdf-corr.eps'

set style line 10 lc 0 lt 0 lw 5
set style line 1  lc 1 lt 1 lw 5 
set style line 2  lc 2 lt 1 lw 5 

set xtics 0.1
set mxtics 10
set ytics .5
set mytics 5

set xlabel "r [nm]"

pl [0:.5][0:2.5] 1 ls 10 not,\
'rdf.ho.xvg' u 1:2 ls 1 w l t 'H-O RDF',\
'rdf.hh.xvg' u 1:2 ls 2 w l t 'H-H RDF'
