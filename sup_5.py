# This code was used for processing nnLuz thermostability data (Supplementary Figure 5) 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D
import scipy.stats as ss
import scikit_posthocs as scihoc
import seaborn as sns
import string

from os import listdir
from re import search

path = '' # Insert here path to raw data
experiment_to_analyse = '' # ID of the experiment

annotations = pd.read_csv(f'', delimiter=';') # Insert here path to file with roi-to-names mapping

anno_dict = {}
for index, record in annotations.iterrows():
  anno_dict[f'Mean{record["Order"]}'] = f'{record["Name"]}-{record["ppas"]}-{record["Spot"]}'


data = pd.read_csv(f'{path}.csv', index_col=' ') # Insert here path to raw data
data_cols = [x for x in data.columns if 'Mean' in x]
data = data[data_cols]
data = data.rename(columns=anno_dict)


not_bg_cols = [x for x in data.columns if 'BG' not in x]
bg_cols = [x for x in data.columns if 'BG' in x]
bg_mean = data[bg_cols].mean(axis = 1)
bg_std = data[bg_cols].std(axis = 1)
bg = bg_mean - 3 * bg_std
for col in data.columns:
  data[col] -= bg
data = data[not_bg_cols]

data = data.transpose()
data['Label'] = data.apply(lambda row: row.name.split('-')[0], axis=1)
data['ppas'] = data.apply(lambda row: row.name.split('-')[1], axis=1)
data['Sum'] = data.apply(lambda row: row[0:30].sum(), axis=1)

data.to_csv(f'{path}{experiment_to_analyse}_data.csv')

# Normalization by values at 25 C
data = pd.read_csv(f'{path}{experiment_to_analyse}_data.csv', usecols = ['Unnamed: 0', 'Label', 'ppas', 'Sum'])
data['Temp'] = data.apply(lambda row: row.Label.split(', ')[1], axis = 1)
for ppas in data.ppas.unique():
  selected = data[data.ppas == ppas]
  norm_value = selected[selected.Temp == '25C']['Sum'].mean()
  for name in selected['Unnamed: 0'].unique():
    index_ = data[data['Unnamed: 0'] == name].index.values[0]
    data.at[index_, 'Sum'] /= norm_value

fig, axes = plt.subplots(1, figsize = (12,6))

sorting_order = []
for luz in ['nnLuz WT', 'nnLuz v4']:
  for temp in ['25C', '30C', '35C', '40C', '45C', '50C']:
    sorting_order.append(f'{luz}, {temp}')

medianprops = {'color': 'coral'}
capprops = {'color': 'black'}
whiskerprops = {'color': 'black'}
boxprops = {'edgecolor': 'black'}

ax1 = axes
sns.boxplot(data=data,
            x='Label',
            y='Sum',
            color = 'white',
            order = sorting_order,
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            showfliers = False,
            ax = ax1)

sns.swarmplot(data=data, x='Label', y='Sum', hue = 'ppas', order = sorting_order, ax = ax1)

ax1.set_yscale('log')
ax1.set_xlabel('Temperature, $^\circ$C', fontsize = 20)
ax1.set_ylabel('Normalized integral\n luminescence, RLU', fontsize = 20)



ax1.grid(lw = 0.5)
ax1.legend().set_visible(False)
sns.despine(offset = 10, trim = False, ax = ax1)
ax1.set_xticks(ticks = np.arange(0, len(sorting_order), 1))
ax1.set_xticklabels(labels = [x.replace('C', '').split(', ')[1] for x in sorting_order])
ax1.tick_params(axis = 'y', labelsize = 15)
ax1.tick_params(axis = 'x', labelsize = 17, rotation = 0)


ax2 = plt.axes([0,0,1,1], facecolor=(0,0,0,0))
ax2.axis('off')
rectangles = {'nnLuz WT': patches.Rectangle((0.113, 0.880), 0.400, 0.063, color = 'whitesmoke'),
              'nnLuz v4': patches.Rectangle((0.513, 0.880), 0.385, 0.063, color = 'silver'),
              }

for r in rectangles:

    ax2.add_patch(rectangles[r])
    rx, ry = rectangles[r].get_xy()
    cx = rx + np.divide(rectangles[r].get_width(), 2.0)
    cy = ry + np.divide(rectangles[r].get_height(), 2.0)

    ax2.annotate(r,
                 (cx, cy),
                 color='black',
                 fontsize=20,
                 ha='center',
                 va='center')


plt.subplots_adjust(hspace=3)
ax1.set_title(r'$\it{Pichia}$' + ' ' + r'$\it{pastoris}$' + ', stable, after 50 \u03BCM luciferin treatment \n\n', fontsize = 22)

plt.savefig(f'', 
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='white',
            pad_inches=0)
