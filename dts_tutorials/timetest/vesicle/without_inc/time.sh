
#echo $FILES
#read -p "If you want to run Beads snapshot input 1 and if you want to run triangulation snapshot input 2: " name

#!/bin/bash

# Line width of the axes
set border linewidth 1 


# Line styles
set style line 1 linecolor rgb 'cyan'   lt 2 linewidth 3 pt 1
set style line 2 linecolor rgb 'red' lt 2 linewidth 3 pt 2
set style line 3 linecolor rgb 'black' lt 2 linewidth 3 pt 3
set style line 4 linecolor rgb 'blue' lt 2 linewidth 3 pt 4
set style line 5 linecolor rgb 'green' lt 2 linewidth 3 pt 5
set style line 6 linecolor rgb 'purple' lt 2 linewidth 3 pt 6

set terminal postscript enhanced color solid eps 16 size 5 ,7  font "Times-Roman, 16" #font 'UTF-8'#font 'Times-Roman' #30 #font 'Verdana,10'

set style line 12 lc rgb '#ddccdd' lt 1 lw 2.5 # --- red


set output "vesicle_L_T_INC.eps"
#set key at 9 , 3
set xlabel 'move [#]'
set ylabel 'time(ms per 10^3 moves)'
# Axes ranges
#set xrange [-5:6] reverse
#set yrange [-50:5]
# Axes tics
#set xtics font "Times-Roman, 30"
#set xtics ('-2π' -2*pi, '-π' -pi, 0, 'π' pi, '2π' 2*pi)
#set ytics 5
#set xtics 2
#set tics scale 0.75
#set mxtics 1
#set mytics 0.1

#set grid xtics mxtics ytics mytics back ls 12, ls 12
set multiplot layout 2,1

plot "link_time_ms.xvg"  u 1:(2*($3-0.5)*$2) t 'link flip'
plot "vertex_time_ms.xvg" u  1:(2*($3-0.5)*$2) t 'vertex move'





