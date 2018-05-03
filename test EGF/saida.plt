set encoding utf8
set title 'EGF concentration -  Agents: 4    Time step: 0'
set dgrid 35,35
set palette rgbformulae -31,13,22
set grid
#set cbrange[0:0.006]
set xtics 0,34,340
set ytics 0,34,340
set view map
set pm3d interpolate 0,0
set term postscript eps enhanced color
set size ratio -1
set output 'egf-00000.eps'
splot 'egf-00000.dat' u 1:2:3 with pm3d notitle
