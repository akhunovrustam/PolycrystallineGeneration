reset
set grid
set grid back
set terminal postscript eps enhanced color font 'Helvetica,10'
#set terminal windows
set output "comparison.eps"
file_exists(file) = system("if exist ".file." echo 1")
set yrange [0:0.37]
set xrange [0:20000]
set xlabel 'Number of evoluation of fitness penalty function' font ',15'
set ylabel 'Value of penalty' font ',15'
set xtics font ", 14"
set ytics font ", 14"
set key font ", 14"

set multiplot layout 1,1 rowsfirst  
set print 'test.txt'
cp = "1p 2p uniform"
md = "uniform normal reinit"
mp = "0.01"
ps = "70"
fac = "1.0"
from = "10"
to = "100"
stat = "best"
step = 1
set object 2 rect from 12500,0.19 to 20000,0.35 fc rgb "#eeeeee"
colors = "black dark-grey red web-green web-blue dark-magenta dark-cyan dark-orange dark-yellow royalblue goldenrod"
do for [i=1:words(cp)] {
	do for [j=1:words(md)] {
		do for [k=1:words(mp)] {
			do for [m=1:words(ps)] {
				do for [fc=1:words(fac)] {
					do for [fr=1:words(from)] {
						do for [t=1:words(to)] {
							do for [st=1:words(stat)] {
								config = word(cp, i).'_'.word(md, j).'_'.word(mp, k).'_'.word(ps, m).'_'.word(fac, fc).'_'.word(from, fr).'_'.word(to, t).'_'.word(stat, st)
								config_l = word(cp, i).'\_'.word(md, j).'\_'.word(mp, k).'\_'.word(ps, m).'\_...'
								#print config
								bla = '2,0.2,0'
								set key inside
								set key right top
								set border 3
								set key at 20000,(0.35-step*0.015)
								fl = config.".csv"
								if  ( file_exists(fl) eq '1' ) { 
									print fl
									plot fl u 2:3:4 w yerrorbars title config_l linetype rgb word(colors, step) lw 2; 
									#plot fl u 2:3 smooth cspline title config_l linetype rgb word(colors, step) lw 2
									print '#00'.i.i.j.j
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



unset multiplot 