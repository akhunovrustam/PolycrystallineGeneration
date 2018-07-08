import glob
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

folder = 'test/results_orientation'
cp = ['1p']
# cp = ['1p', '2p']
# md = ['uniform', 'normal', 'reinit']
md = ['uniform']
mp = ['0.8']
# mp = ['0.8', '0.4', '0.2', '0.1', '0.05', '0.02', '0.01']
# ps = ['10', '30', '50', '70', '90', '110']
ps = ['10', '90']
# mult = ['1.0', '2.0', '4.0', '8.0', '16.0']
mult = ['1.0']
# fro = ['10', '20']
fro = ['10']
to = ['200']
# to = ['200']
# thr1 = ['50', '100', '200', '300', '400', '500']
thr1 = ['200']
thr2 = ['100']
# thr2 = ['100', '200', '300', '500', '600', '700']
fac = ['10']
ps_fac1 = ['1']
ps_fac2 = ['1']

# stats = ['avg', 'best', 'worst']
stats = ['best']
delim = "_"

plt.rcParams["figure.figsize"] = [16,6]
# fig = plt.figure()
ax = plt.subplot(111)
plt.xticks(rotation=45)
plt.xlabel('Number of penalty evaluation')
plt.ylabel('Penalty')
plt.subplots_adjust(bottom=0.16)
Step = 300

for i in cp:
	for j in md:
		for k in mp:
			for m in ps:
				for mu in mult:
					for f in fro:
						for t in to:
							for t1 in thr1:
								for t2 in thr2:
									for fc in fac:
										for pfc1 in ps_fac1:
											for pfc2 in ps_fac2:
												# config = '2p_reinit_0.8_100'
												# config = i + delim + j + delim + k + delim + m + delim + mu + delim + f + delim + t
												config = i + delim + j + delim + k + delim + m + delim + mu + delim + f + delim + t \
													+ delim + t1 + delim + str(int(t1) + int(t2)) + delim + fc + delim + pfc1 + delim + pfc2
												print(config)
												
												filelist = glob.glob(folder + '/' + config + '_SEED*')
												cnt = 0
												min = 1000000
												if len(filelist) == 0: continue
												for file in filelist:
													dt = pd.read_csv(file + '\penalty_steps_avg.txt', delim_whitespace=True)
													if min > len(dt.index): min = len(dt.index)
													
												for st in stats:
													common_frame = pd.DataFrame
													for file in filelist:
													# print(file)
														try:
															dt = pd.read_csv(file + '\penalty_steps_' + st + '.txt', delim_whitespace=True)
															if min > len(dt.index): min = len(dt.index)
															if common_frame.empty:
																x = np.append([0], np.array(dt['Step'][int(100/dt['Step'][0] - 1):min:int(Step/dt['Step'][0])]))
																y = x.astype(np.str)
																common_frame = pd.DataFrame(columns=y);
																print(file)
																common_frame.loc[cnt] = np.append([dt['Penalty'][0]], 
																	np.array(dt['Penalty'][int(100/dt['Step'][0] - 1):min:int(Step/dt['Step'][0])]))
																cnt += 1;
															else:
																print(file)
																common_frame.loc[cnt] = np.append([dt['Penalty'][0]], 
																	np.array(dt['Penalty'][int(100/dt['Step'][0] - 1):min:int(Step/dt['Step'][0])]))
																cnt += 1;
														except TypeError:
															continue
													result_frame = pd.DataFrame(columns=['Step', 'Avg Penalty', 'Std deviation'])
													
													# if config == '2p_reinit_0.8_100': print(common_frame[col], common_frame[col].map(type).unique(), len(common_frame[col]))
													
													for num in range(0, len(common_frame.columns)):
														col = ''
														col = common_frame.columns[num]
														result_frame.loc[num] = [col, common_frame[col].mean(), common_frame[col].std()]
														
													try:
														# print(result_frame);
														result_frame.to_csv(folder + '\\' + config + '_' + st + '.csv', sep=' ')
														# if st == 'avg':
															# ax.plot('Step', 'Avg Penalty', data=result_frame, marker='o', label="")
														# elif st == 'best':
															# ax.plot('Step', 'Avg Penalty', data=result_frame, marker='^', label="")
														# elif st == 'worst':
															# ax.plot('Step', 'Avg Penalty', data=result_frame, marker='v', label="")
														ax.errorbar('Step', 'Avg Penalty', 'Std deviation', data=result_frame, marker='.', label=config + " " +st)
													# break;
													except ValueError:
														continue
												# plt.yscale('log')
												# ax.title('Step - Penalty')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=5, fancybox=True, shadow=False, framealpha=0.5)
plt.legend()
plt.show()
				# exit()