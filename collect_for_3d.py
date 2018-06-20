import glob
import pandas as pd
import math
import numpy as np
import numbers

folder = 'test/results_size'
# cp = ['2p']
cp = ['1p', '2p', 'uniform']
# md = ['uniform', 'normal', 'reinit']
md = ['uniform']
mp = ['0.8', '0.2', '0.1', '0.05', '0.02', '0.01']
# mp = ['0.1']
ps = ['10', '30', '50', '70', '90', '110']
# ps = ['10']
# mult = ['1.0', '2.0', '4.0', '8.0', '16.0']
mult = ['1.0']
# fro = ['10', '20']
fro = ['10']
# to = ['100']
to = ['100']

# stats = ['avg', 'best', 'worst']
stats = ['best']
delim = "_"

rframe = pd.DataFrame(columns=['mp', 'ps', 'penalty'])
num = 0
for i in cp:
	for j in md:
		for k in mp:
			for m in ps:
				for mu in mult:
					for f in fro:
						for t in to:
							# config = '2p_reinit_0.8_100'
							config = i + delim + j + delim + k + delim + m + delim + mu + delim + f + delim + t
							filelist = glob.glob(folder + '/' + config + '_SEED*')
							cnt = 0
							min = 1000000
							if len(filelist) == 0: continue
							print(config)
							avg = 0
							cnt = 0
							for file in filelist:
								dt = pd.read_csv(file + '\penalty_steps_best.txt', delim_whitespace=True)
								if isinstance(dt['Penalty'].iloc[-1], numbers.Number):
									avg += dt['Penalty'].iloc[-1]
									cnt += 1
									
							avg = avg / cnt;
							rframe.loc[num] = [k, m, avg]
							num += 1
							print(avg)
rframe.to_csv('mp_ps_iso.csv', sep=' ')