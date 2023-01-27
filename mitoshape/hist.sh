
#echo $FILES
#read -p "If you want to run Beads snapshot input 1 and if you want to run triangulation snapshot input 2: " name

#!/bin/bash

# Line width of the axes
set border linewidth 2 


# Line styles
set style line 1 linecolor rgb 'red'   lt 1 linewidth 2 pt 7 ps 2.5
set style line 2 linecolor rgb 'blue' lt 1 linewidth 2 pt 1 ps 5
set style line 3 linecolor rgb 'blue' lt 1 linewidth 5 pt 6 ps 2.5

set style line 4 linecolor rgb 'blue'   lt 1 linewidth 5 pt 1
set style line 5 linecolor rgb 'purple' lt 1 linewidth 5 pt 0
set style line 6 linecolor rgb 'purple' lt 1 linewidth 5 pt 1


#set style line 4 linecolor rgb '#00FF00' lt 3 linewidth 1 pt 9
#set style line 5 linecolor rgb '#808000' lt 2 linewidth 1 pt 13
#set style line 6 linecolor rgb '#A52A2A' lt 3 linewidth 1 pt 9
#set style line 7 linecolor rgb '#FFD700' lt 2 linewidth 1 pt 13
#set style line 8 linecolor rgb '#008080' lt 3 linewidth 1 pt 9


#set terminal epslatex size 9cm,7cm color colortext standalone header \
#"\\newcommand{\\ft}[0]{\\footnotesize}"


#set terminal postscript enhanced color solid eps 16 size 4 , 3 font 'Times-Roman' #30 #font 'Verdana,10'


#set terminal postscript postscript      # old postscript
#set terminal postscript enhanced      # old enhpost
#set terminal postscript landscape 22  # old psbig
#set terminal postscript eps 14        # old epsf1
#set terminal postscript eps 14       # old epsf2
#set size 0.7,1.4; set term post portrait color "Times-Roman" 14


set terminal postscript enhanced color solid eps 14 size 6 ,3 #font 'UTF-8'#font 'Times-Roman' #30 #font 'Verdana,10'





#set label 2 '\ft $5$\,meV'        # at 1     rotate by  78.5 center tc ls 1
#set label 3 '\ft $10$\,meV'       # at 2   rotate by  71.8 center tc ls 2

set style line 12 lc rgb '#ddccdd' lt 1 lw 2.5 # --- red


set output "curvehist.eps"
#set key at 9 , 3
set ylabel 'Sensing'
set xlabel 'C_{0}R_{tube}'
# Axes ranges
#set xrange [0:1]
#set yrange [-1:5]
# Axes tics
#set xtics font "Times-Roman, 30"
#set xtics ('-2π' -2*pi, '-π' -pi, 0, 'π' pi, '2π' 2*pi)
set ytics 0.05
#set xtics 0.2
#set ytics scale 0.01set ylabel 'Density [#/vertex]'
#set mxtics 0.01
#set mytics 0.1

#set grid xtics mxtics ytics mytics back ls 12, ls 12
set multiplot layout 1,2





######## Plot 1 DOPC sn1 order parameter
#unset key
#unset xlabel
#unset ylabel


######## Plot 2

set key  #at -0.2,5
#set key font #", 15"
#set xlabel  '{Carbon number}'
#unset ylabel
set title ""
#set label 1 '(B)' at graph 0.05,0.92 #font 'Times-Roman'
#set lmargin screen 0.2
#set rmargin screen 0.95
#set bmargin at screen 0.2
#set tmargin at screen 0.9
#unset ytics'
#set logscale y
#set logscale x
#unset key
#unset ytics
#set key default
#set key box
#set key box linestyle 1
#set key out vert
#:3 with yerrorbars

#set lmargin screen 0.3
#set rmargin screen 0.9

#set label 1 '{/Symbol t} = 1' at graph 0.05,0.9 font ',20'
#set ylabel '{/Symbol k}_{eff}[k_BT]' #offset 0,0.4,0
#set xlabel 'coverage ' #offset 1.7,0,0

#unset key

#unset key

#set xrange [0:1.1]
#set yrange [0:0.44]
#set lmargin screen 0.25
#set rmargin screen 0.95
#set bmargin at screen 0.2
#set tmargin at screen 0.9
#unset key
#set label 1 '{/Symbol D}{/Symbol k} = 15, {/Symbol D}{/Symbol k}_{G} = 10' at graph 0.05,0.9 font ',20'

#plot "data.xvg" u ($1*sqrt(5/2)):2  ls 1 t '{/Symbol D}{/Symbol k} = 15, {/Symbol D}{/Symbol k}_{G} = 10' smooth sbezier,\
#"data.xvg" u ($1*sqrt(5/2)):3  ls 3 t '{/Symbol D}{/Symbol k} = 30, {/Symbol D}{/Symbol k}_{G} = 25' smooth sbezier ,\

#unset xtic
#unset ytic
#unset ylabel
#unset key
set key box width 2 height 2 opaque
set xrange [-0.2:0.2]
plot "histogram_mean_curvature.xvg" u 1:2  w l ls 1 t "mean"
set xrange [-0.02:0.02]
plot "histogram_Gauss_curvature.xvg" u 1:2 w l  ls 2 t "gauss"


#, "R8-tether.xvg" u (11.2*$1/1000):($4)  ls 4 t "mixed: C_0 = 8.7  {/Symbol m} m^{-1} 50 25 ",
#"R9-tether.xvg" u (11.2*$1/1000):($4)  ls 6 t "mixed: C_0 = 8.7  {/Symbol m} m^{-1} 80 50 "

#unset ytic
#unset ylabel
#set label 1 '{/Symbol t} = 1' at graph 0.05,0.9 font ',20'
unset key

set ytic

unset key
unset label 1
#set ylabel 'R [d]'
unset xlabel
unset xtic
#set yrange [0:0.09]

#set ylabel 'R(z)'


#unset ylabel
#unset ytic
set key # at 170,2.5 font ',13'
set key  font ',15'







