set multiplot layout 2,1 rowsfirst                                 
plot "dist_first.txt" u 1:2 w boxes, "dist_first.txt" u 1:3 w lines
plot "dist_tmp.txt" u 1:2 w boxes, "dist_first.txt" u 1:3 w lines  
unset multiplot 