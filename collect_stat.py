import glob
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

folder = 'results'
cp = ['1p']
# cp = ['1p', '2p', 'uniform']
# md = ['uniform', 'normal', 'reinit']
md = ['uniform']
# mp = ['0.8' '0.2' '0.05' '0.01']
mp = ['0.8']
ps = ['10', '50', '100']

delim = "_"

plt.rcParams["figure.figsize"] = [12,5]
# fig = plt.figure()
ax = plt.subplot(111)
plt.xticks(rotation=90)

for i in cp:
	for j in md:
		for k in mp:
			for m in ps:
				# config = '2p_reinit_0.8_100'
				config = i + delim + j + delim + k + delim + m
				filelist = glob.glob(folder + '/' + config + '_SEED*')
				common_frame = pd.DataFrame
				cnt = 0
				min = 1000000
				if len(filelist) == 0: continue
				print(config)
				for file in filelist:
					dt = pd.read_csv(file + '\penalty_steps_avg.txt', delim_whitespace=True)
					if min > len(dt.index): min = len(dt.index)
					
				for file in filelist:
					# print(file)
					dt = pd.read_csv(file + '\penalty_steps_avg.txt', delim_whitespace=True)
					if min > len(dt.index): min = len(dt.index)
					if common_frame.empty:
						x = np.array(dt['Step'][0:int(min/32)])
						y = x.astype(np.str)
						common_frame = pd.DataFrame(columns=y);
						print(file)
						common_frame.loc[cnt] = np.array(dt['Penalty'][0:int(min/32)])
						cnt += 1;
					else:
						print(file)
						common_frame.loc[cnt] = np.array(dt['Penalty'][0:int(min/32)])
						cnt += 1;
						
				result_frame = pd.DataFrame(columns=['Step', 'Avg Penalty', 'Std deviation'])
				
				# if config == '2p_reinit_0.8_100': print(common_frame[col], common_frame[col].map(type).unique(), len(common_frame[col]))
				
				for num in range(0, len(common_frame.columns)):
					col = ''
					col = common_frame.columns[num]
					result_frame.loc[num] = [col, common_frame[col].mean(), common_frame[col].std()]
					
				try:
					# print(result_frame);
					result_frame.to_csv(folder + '\\' + config + '.csv', sep=' ')
					ax.plot('Step', 'Avg Penalty', data=result_frame, marker='x', label="")
					ax.errorbar('Step', 'Avg Penalty', 'Std deviation', data=result_frame, marker='.', label=config)
				# break;
				except ValueError:
					continue
				# plt.yscale('log')
				# ax.title('Step - Penalty')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=5, fancybox=True, shadow=False, framealpha=0.5)
# plt.legend()
plt.show()
				# exit()