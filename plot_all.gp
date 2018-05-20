file_exists(file) = system("if exist ".file." echo 1")

set multiplot layout 1,1 rowsfirst  
set print 'test.txt'
cp = "1p 2p uniform"
md = "uniform normal reinit"
mp = "0.8"
ps = "10 50 100"
do for [i=1:words(cp)] {
	do for [j=1:words(md)] {
		do for [k=1:words(mp)] {
			do for [m=1:words(ps)] {
				config = word(cp, i).'_'.word(md, j).'_'.word(mp, k).'_'.word(ps, m)
				config_l = word(cp, i).'\_'.word(md, j).'\_'.word(mp, k).'\_'.word(ps, m)
				#print config
				bla = '2,0.2,0'
				set key at -1.0,1.0
				fl = config.".csv"
				if  ( file_exists(fl) eq '1' ) { 
					print fl
					plot fl u 2:3:4 w yerrorbars title config_l linetype i+j+k+m; 
				}
			}
		}
	}
}
unset multiplot 