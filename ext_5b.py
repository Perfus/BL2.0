# This code was used for processing and plotting data for Extended Data 5b

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import scipy.stats as ss
import scikit_posthocs as scihoc
import string

from matplotlib import patches
from matplotlib.lines import Line2D
from re import search

path = '' # Insert here path to raw data
exp = '' # Insert here experiment ID

filelist = [x.replace('_roi_names.xlsx', '') for x in os.listdir(path) if 'roi_names' in x]

def stat_rounder(p):
  if p < 0.0001:
    return 'p < 0.0001'
  else:
    return f'p={str(round(p, 4))}'

def data_prepare(filename, savename_suffix):
  data = pd.read_csv(f'{path}{filename}', index_col = ' ')

  data['long_label'] = data.apply(lambda row: row.name, axis = 1)
  data['long_label'] = data.apply(lambda row: anno_dict[row['long_label']], axis = 1)
  data.set_index('long_label', inplace = True)
  bg_rows = [x for x in data.index if 'BG' in x]
  bg = data[data.index.isin(bg_rows)].mean()

  for index, row in data.iterrows():
    data.loc[index] -= bg

  data = data[~data.index.isin(bg_rows)]
  data['Label'] = data.apply(lambda row: row.name.split('__')[0], axis = 1)
  data['Plant'] = data.apply(lambda row: row.name.split('__')[1], axis = 1)
  data.to_csv(f'{path}{exp_name}_{savename_suffix}_data_transformed.csv')

  return data

data = pd.DataFrame()
for file in filelist:
  annotations = pd.read_excel(f'{path}{file}_roi_names.xlsx')
  anno_dict = dict((row['Order'], f"{row['Plasmid']}__{row['Plant']}__{row['Order']}") for index, row in annotations.iterrows())

  prep = data_prepare(f'{file}.csv', file)

  data = pd.concat([data, prep])

data.to_csv(f'{path}{exp_name}full_data.csv')
mann_whitney_test = scihoc.posthoc_mannwhitney(data, val_col = 'Mean', group_col = 'Label', p_adjust = 'holm-sidak')

medianprops = {'color': 'coral'}
capprops = {'color': 'black'}
whiskerprops = {'color': 'black'}
boxprops = {'edgecolor': 'black'}

xlabel_size = 20
ylabel_size = 20
yticklabel_size = 15
xticklabel_size = 14
title_size = 24
suptitle_size = 20
legend_size = 20
signplot_size = 16
letter_size = 24

fig, axes = plt.subplots(1, figsize = (4, 6))


sorting_order = data.groupby('Label').median().sort_values('Mean', ascending = True).index.values

ax1 = axes
sns.boxplot(data=data,
            x='Label',
            y='Mean',
            color = 'white',
            order = sorting_order,
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            ax = ax1)

sns.swarmplot(data=data, x='Label', y='Mean', hue='Plant', order = sorting_order, ax = ax1)

ax1.set_yscale('log')
ax1.set_xlabel(None)
ax1.set_ylabel('Luminescence, RLU', fontsize = ylabel_size)
ax1.set_ylim(100, 1100)

ax1.grid(lw = 0.5)
ax1.legend().set_visible(False)
sns.despine(offset = 10, trim = False, ax = ax1)
ax1.set_xticks(ticks = np.arange(0, len(sorting_order), 1))
ax1.tick_params(axis = 'x', labelsize = xticklabel_size)
ax1.set_xticklabels(labels = [x for x in sorting_order])
ax1.tick_params(axis = 'y', labelsize = yticklabel_size)

y1 = 150
# statistical annotation
x1, x2 = 0, 1
y, h, col = y1, y1/10, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data[data["Label"] == "FBP3"]["Mean"].mean(), data[data["Label"] == "FBP2"]["Mean"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test["FBP3"]["FBP2"])}'
ax1.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

plt.suptitle(r'$\it{Nicotiana}$' + ' ' + r'$\it{benthamiana}$' + ',\nstable lines', fontsize = 22)

plt.savefig('', # Insert here path for figure to save
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='white')
