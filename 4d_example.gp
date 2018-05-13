set key off
set xrange [-1:2]
set yrange [-1:2]
set zrange [0:0.2]
set multiplot layout 2,1 rowsfirst                                 
splot 'dist_first.txt' using 1:2:3:4 with points palette pointsize 1 pointtype 7
splot 'dist_tmp.txt' using 1:2:3:4 with points palette pointsize 1 pointtype 7
unset multiplot 