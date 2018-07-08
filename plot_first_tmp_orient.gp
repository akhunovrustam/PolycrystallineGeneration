reset
set xrange [0:70]
set yrange [0:0.05]
#set terminal postscript eps enhanced color font 'Helvetica,10'
#set output 'first_last_comparison.eps'
set multiplot layout 2,1 rowsfirst
set key inside
set ylabel 'Normalized frequency'
set xlabel 'Misorientation angle'
plot "dist_first.txt" u 1:2 w boxes title 'Random orientation', "dist_tmp.txt" u 1:3 w lines title 'Mackenzie'
plot "dist_tmp.txt" u 1:2 w boxes title 'After calculation', "dist_tmp.txt" u 1:3 w lines title 'Mackenzie'
unset multiplot 