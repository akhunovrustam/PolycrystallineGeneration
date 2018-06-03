reset
set xrange [0:2.5]
set yrange [0:3]
#set terminal postscript eps enhanced color font 'Helvetica,10'
#set output 'first_last_comparison.eps'
set multiplot layout 2,1 rowsfirst
set key inside
set ylabel 'Normalized frequency'
set xlabel 'Normalized grain size (D/<D>)'
plot "dist_first.txt" u 1:2 w boxes title 'Random tessellation', "dist_first.txt" u 1:3 w lines title 'Lognormal distribution'
plot "dist_tmp.txt" u 1:2 w boxes title 'Tessellation after calculations', "dist_first.txt" u 1:3 w lines title 'Lognormal distribution', \
	"dist_tmp.txt" u 1:4 w points title 'Error', "dist_tmp.txt" u 1:5 w lines title 'Error sum
unset multiplot 