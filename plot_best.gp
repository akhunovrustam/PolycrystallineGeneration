reset
set grid
set grid back
#set terminal postscript eps enhanced color font 'Helvetica,10'
#set output "comparison.eps"
set terminal windows
set yrange [0:0.4]
set xrange [0:45500]
set xlabel 'Number of evoluation of fitness penalty function' font ',15'
set ylabel 'Value of penalty' font ',15'
set xtics font ", 14"
set ytics font ", 14"
set key font ", 14"

set print 'test.txt'
config = '1p_uniform_0.01_70_1.0_10_200_200_500_10_1_1'
config_l = 'The best configuration'

print config
set key inside
set key right top
set border 3
set key at 45000,(0.40-step*0.02)
#plot fl u 2:3:4 w yerrorbars title config_l linetype rgb word(colors, step) lw 2; 
plot config.'_avg.csv' u 2:3 w l title config_l.' average result' lw 2, config.'_best.csv' u 2:3 w l title.' best result' config_l lw 2, \
	config.'_worst.csv' u 2:3 w l title.' worst result' config_l lw 2

set terminal windows