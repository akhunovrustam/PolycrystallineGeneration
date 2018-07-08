reset
set grid
set grid back
set terminal postscript eps enhanced color font 'Helvetica,10'
#set terminal windows
set output "comparison.eps"
file_exists(file) = system("if exist ".file." echo 1")
set yrange [0:0.37]
set xrange [0:20500]
set xlabel 'Number of evoluation of fitness penalty function' font ',15'
set ylabel 'Value of penalty' font ',15'
set xtics font ", 14"
set ytics font ", 14"
set key font ", 14"

set multiplot layout 1,1 rowsfirst  
set print 'test.txt'
#cp = "1p 2p uniform"
cp = "1p"
#md = "uniform normal reinit"
md = "uniform"
mp = "0.01"
#ps = "10 30 50 70 90 110"
ps = "70"
fac = "1.0"
from = "10"
to = "200"
#thr1 = "50 100 200 300 400 500"
thr1 = "200"
thr2 = "100 200 300 500 600 700"
#thr2 = "100"
mpfac = "10"
thr1fac = "1"
thr2fac = "1"
stat = "best"
step = 1
set object 2 rect from 13200,0.25 to 20000,0.35 fc rgb "#eeeeee"
colors = "black dark-grey red web-green web-blue dark-magenta dark-cyan dark-orange dark-yellow royalblue goldenrod"
do for [i=1:words(cp)] {
	do for [j=1:words(md)] {
		do for [k=1:words(mp)] {
			do for [m=1:words(ps)] {
				do for [fc=1:words(fac)] {
					do for [fr=1:words(from)] {
						do for [t=1:words(to)] {
							do for [th1=1:words(thr1)] {
								do for [th2=1:words(thr2)] {
									do for [mfc=1:words(mpfac)] {
										do for [th1f=1:words(thr1fac)] {
											do for [th2f=1:words(thr2fac)] {
												do for [st=1:words(stat)] {
													config = word(cp, i).'_'.word(md, j).'_'.word(mp, k).'_'.word(ps, m).'_'.word(fac, fc).\
														'_'.word(from, fr).'_'.word(to, t).\
														'_'.word(thr1, th1).'_'.(word(thr1, th1)+word(thr2, th2)).'_'.word(mpfac, mfc).\
														'_'.word(thr1fac, th1f).'_'.word(thr2fac, th2f).'_'.word(stat, st)
													config_l = '...\_'.word(from, fr).'\_'.word(to, t).\
														'\_'.word(thr1, th1).'\_'.(word(thr1, th1)+word(thr2, th2)).'\_...'
													print config
													bla = '2,0.2,0'
													set key inside
													set key right top
													set border 3
													set key at 20000,(0.36-step*0.015)
													fl = config.".csv"
													if  ( file_exists(fl) eq '1' ) { 
														print fl
														#plot fl u 2:3:4 w yerrorbars title config_l linetype rgb word(colors, step) lw 2; 
														plot fl u 2:3 w l title config_l linetype rgb word(colors, step) lw 2
													}
													
													if (step == 1) {
														unset object 2
														unset grid
													}
													
													step = step + 1
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}



unset multiplot 