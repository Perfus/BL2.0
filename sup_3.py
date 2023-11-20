# This code was used for processing thermostability in E.coli data (Supplementary Figure 3)

import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import scipy.stats as ss
import scikit_posthocs as scihoc
import seaborn as sns

path = '' # Insert here path to raw data

data = pd.read_excel(f'{path}', index_col = 'time(min)') # Insert here path to raw data
data.rename(columns = {'nLuz WT': 'nnLuz WT_1',
                       'Unnamed: 2': 'nnLuz WT_2',
                       'nLuz Triple': 'nnLuz v3_1',
                       'Unnamed: 4': 'nnLuz v3_2'}, inplace = True)

xlabel_size = 20
ylabel_size = 20
yticklabel_size = 15
xticklabel_size = 16
title_size = 24
suptitle_size = 20
legend_size = 20
signplot_size = 16
letter_size = 24

colors_ = ['darkblue',
          'darkred']

def model_func(t, A, K, C):
    return (A - C) * np.exp(-K * t) + C

def fit_exp_nonlinear(t, y):
    opt_parms, parm_cov = sp.optimize.curve_fit(model_func, t, y)
    A, K, C = opt_parms
    return A, K, C

def model_func(t, A, K, C):
  return (A - C)*np.exp(-K*t) + C

data = data.transpose()
data['Label'] = data.apply(lambda row: row.name.split('_')[0], axis = 1)
fig, axes = plt.subplots(1, figsize = (10, 6))
ax = axes
sorting_order = ['nnLuz WT', 'nnLuz v3']

to_plot = data
time = np.array(data.columns[:-1].values.tolist())


for index, label in enumerate(sorting_order):
  mean = to_plot[to_plot.Label == label].iloc[:, 0:-1].mean(axis = 0).values.tolist()
  std = to_plot[to_plot.Label == label].iloc[:, 0:-1].std(axis = 0).values.tolist()
  A, K, C = fit_exp_nonlinear(time, mean)
  fit_y = model_func(time, A, K, C)
  ax.plot(time, fit_y, color = colors_[index], lw = 3)

  half_time = np.log(2)/K

  plot = ax.plot(time, to_plot[to_plot.Label == label].iloc[0, 0:-1], color = colors_[index], label = f'{label}, T' + u'\u2080.\u2085' + f'= {round(half_time, 1)} min', marker = 's',
                   lw = 0, markersize = 7, mec = 'black' )
  # plot[-1][0].set_linestyle('--')

  plot = ax.plot(time, to_plot[to_plot.Label == label].iloc[1, 0:-1], color = colors_[index], marker = 's',
                   lw = 0, markersize = 7, mec = 'black' )
  # plot[-1][0].set_linestyle('--')
  # plot = ax.errorbar(x = time, y = mean, yerr = std, color = colors_[index], label = f'{label}, T' + u'\u2080.\u2085' + f'= {round(half_time, 1)} min', marker = 's',
  #                    lw = 2,fmt='.', markersize = 7 , capsize = 9)


ax.legend(loc = 1, fontsize = legend_size, frameon = False, )
ax.set_xlabel('Time, min', fontsize = xlabel_size)
ax.set_ylabel('Normalized luminescence, RLU', fontsize = ylabel_size)

ax.set_xlim(0)
ax.set_ylim(0)
ax.grid(lw = 0.4)
ax.tick_params(axis = 'both', labelsize = xticklabel_size)


plt.subplots_adjust(hspace = 1, wspace=1.2)
plt.suptitle(r'$\it{Escherichia}$' + ' ' + r'$\it{coli}$' + f', lysates,\nafter 50 \u03BCM luciferin treatment', fontsize = title_size)
sns.despine(offset = 10, trim = False, ax = ax)

plt.savefig(f'', # Insert here path for figure to save
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='white')
