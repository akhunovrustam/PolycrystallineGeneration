import pandas as pd
import matplotlib.pyplot as plt

bla1 = pd.read_csv('c:/Users/akhun/Desktop/diploma/results/2p_normal_0.1_4_exp11-5-2018_11.48/penalty_steps_avg.txt', delim_whitespace=True)
# bla.plot(x='Step', y='Penalty')
# plt.title('AVG')

bla2 = pd.read_csv('c:/Users/akhun/Desktop/diploma/results/2p_normal_0.1_4_exp11-5-2018_11.48/penalty_steps_best.txt', delim_whitespace=True)
# bla.plot(x='Step', y='Penalty')
# plt.title('Best')

bla3 = pd.read_csv('c:/Users/akhun/Desktop/diploma/results/2p_normal_0.1_4_exp11-5-2018_11.48/penalty_steps_worst.txt', delim_whitespace=True)
# bla.plot(x='Step', y='Penalty')
# plt.title('Worst')

plt.plot('Step', 'Penalty', data=bla1, marker='o')
plt.plot('Step', 'Penalty', data=bla2, marker='+')
plt.plot('Step', 'Penalty', data=bla3, marker='x')
plt.yscale('log')
plt.show()