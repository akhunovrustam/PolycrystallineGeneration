reset
set xrange [0:2.5]
set yrange [0:3]
set terminal pngcairo
if (!exists("filename")) filename='first_last_comparison.png'
set output sprintf("%s%s%s", prefix, "/", filename)
set multiplot layout 2,1 rowsfirst
set key inside
set label filename at 0.1,1.8
set label penalty at 0.1,2.5
set ylabel 'Normalized frequency'
set xlabel 'Normalized grain size (D/<D>)'
plot sprintf("%s%s%s", prefix, "/", "dist_first.txt") u 1:2 w boxes title 'Random tessellation', \
	sprintf("%s%s%s", prefix, "/", "dist_first.txt") u 1:3 w lines title 'Lognormal distribution'
plot sprintf("%s%s%s", prefix, "/", "dist_tmp.txt") u 1:2 w boxes title 'Tessellation after calculations', \
	sprintf("%s%s%s", prefix, "/", "dist_first.txt") u 1:3 w lines title 'Lognormal distribution', \
	sprintf("%s%s%s", prefix, "/", "dist_tmp.txt") u 1:4 w points title 'Error'

unset multiplot
