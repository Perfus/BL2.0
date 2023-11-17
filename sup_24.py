# This code was used to plot yeast growth kinetics (Supplementary Fig 24)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as ss
import scikit_posthocs as scihoc
import seaborn as sns

from matplotlib import patches
from matplotlib.lines import Line2D
from re import sub, match, search

colors = ['darkgreen', 
          'deepskyblue',
          'magenta',
          'teal',
          'lightcoral',
          'gray',
          'gold',
          'deepskyblue',
          'slategray',
          'darkred']
FBP = {r'$\bf{WT}$': ['ppasWT'],
       r'$\bf{FBP2}$' + '\nnnHispS\nnnH3H v2\nnnLuz v4\nnnCPH\nNpgA': ['ppas1227', 'ppas1228'],
       r'$\bf{FBP3}$' + '\nmcitHispS\nnnH3H v2\nnnLuz v4\nnnCPH\nNpgA' : ['ppas1219', 'ppas1222']}
def to_fbp(value):
  for key in FBP.keys():
    if value in FBP[key]:
      return key

ygr = pd.read_csv('', index_col = 'Unnamed: 0') # Insert here path to raw data file
ygr['FBP'] = ygr.apply(lambda row: to_fbp(row.Label), axis = 1)

time = [30 * x for x in range(len(ygr.columns[:-7]))] # 

fig = plt.figure(figsize = (12, 9))
for index, label in enumerate(FBP.keys()):
  mean = ygr[ygr.FBP == label].iloc[:, :-7].mean(axis = 0, skipna = True)
  std = ygr[ygr.FBP == label].iloc[:, :-7].std(axis = 0, skipna = True)
  plt.errorbar(time, mean, std, color = colors[index], label = label, errorevery = index + 1)
plt.legend(loc = 'best', fontsize = 18, ncol = 3, frameon = False)
plt.ylabel('OD', fontsize = 20)
plt.xlabel('Time, h', fontsize = 20)
sns.despine(offset = 10, trim = False)
plt.tick_params('both', labelsize = 15)
plt.grid(lw = 0.3)
plt.xlim(0)
ticks = np.arange(0, np.max(time), 180)
plt.xticks(ticks = ticks, labels = [str(int(x/60)) for x in ticks])

plt.savefig('', dpi = 900, bbox_inches = 'tight') # Insert here path for figure to save
