# This code was used for processing data in Supplementary Figure 9b

import string
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
import scipy.stats as ss
import scikit_posthocs as scihoc

path = f'' # Insert path to raw data file

def stat_rounder(p):
  if p < 0.0001:
    return 'p < 0.0001'
  else:
    return f'p={str(round(p, 4))}'

def detect_table(data):
  start = None
  end = None
  for index, val in enumerate(data.iloc[:, 0]):
      if val == 'Cycle Nr.' and data.iloc[index + 1, 1] != 'None':
          start = index
  for index, val in enumerate(data.iloc[start:, 0]):
    if val == 'End Time':
      end = index + start - 1

  if start:
      data_r = data.iloc[start + 1:end, :]
      data_r.columns = data.iloc[start]
      data_r = data_r.set_index('Cycle Nr.')

      return data_r

def detect_od_table(data):
  start = None
  for index, val in enumerate(data.iloc[:, 0]):
    if val == '<>' and data.iloc[index + 1, 1] != 'None':
        start = index
  if start:
      data_r = data.iloc[start + 1:start + 9, :13]
      data_r.columns = data.iloc[start]
      data_r = data_r.set_index('<>')

      return data_r

def plate_process(filename, od_filename, anno_dict, od_control):
  df = pd.read_excel(f'{path}{filename}')
  df = detect_table(df)
  df = df[['Time [s]'] + anno_dict.Order.tolist()]

  df_od = pd.read_excel(f'{path}{od_filename}')
  df_od = detect_od_table(df_od)

  od_data = {}
  for index, row in df_od.iterrows():
    for col in row.index:
      if not math.isnan(row[col]):
          od_data[index+str(col)] = row[col]
  od_control_values = [val for key, val in od_data.items() if key in od_control]
  od_control_values_mean = np.mean(od_control_values)
  od_values = dict((key, np.divide(0.7, val - od_control_values_mean)) for key, val in od_data.items() if key not in od_control)
  print(od_control_values)
  print(od_values)
  print(df.transpose())
  for col in df.columns[1:]:
    df[col] *= od_values[col]

  df = df.transpose()

  df['Label'] = df.apply(lambda row: ''.join(anno_dict[anno_dict.Order == row.name]['Name'].values).replace('_ppas_opt', '').replace('_', ' ').replace('wt', 'WT'), axis = 1)
  df['ppas'] = df.apply(lambda row: ''.join(anno_dict[anno_dict.Order == row.name]['ppas'].values), axis = 1)

  return df

anno_dict = pd.read_csv(f'{path}.csv', delimiter = ';') # Insert here path to roi-to-names mapping file
anno_dict_od = pd.read_csv(f'{path}.csv', delimiter = ';') # Insert here path to roi-to-names mapping file for OD


ppas_data = plate_process('', # Insert here path to raw data file
                          '', # Insert here path to raw data file with OD
                          anno_dict,
                          ['C1', 'C2'])

ppas_data = ppas_data[ppas_data.columns[:50].tolist() + ['Label', 'ppas']]
ppas_data['Sum'] = ppas_data.apply(lambda row: row[:-2].sum(), axis = 1)
ppas_data.to_csv(f'{path}.csv') # Insert here path for processed data to save
ppas_data

medianprops = {'color': 'coral'}
capprops = {'color': 'black'}
whiskerprops = {'color': 'black'}
boxprops = {'edgecolor': 'black'}

xlabel_size = 22
ylabel_size = 22
yticklabel_size = 18
xticklabel_size = 18
title_size = 24
suptitle_size = 22
legend_size = 20
signplot_size = 15

data = ppas_data
mann_whitney_test = scihoc.posthoc_mannwhitney(data, val_col = 'Sum', group_col = 'Label', p_adjust = 'holm-sidak')
sorting_order = ['nnLuz WT', 'nnLuz v4']

fig, axes = plt.subplots(1, figsize = (6, 7))

print(sorting_order)
ax1 = axes
sns.boxplot(data=data,
            x='Label',
            y='Sum',
            color = 'white',
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            order = sorting_order,
            showfliers = False,
            ax = ax1)

sns.swarmplot(data=data, x='Label', y='Sum', hue = 'ppas', order = sorting_order, ax = ax1)

ax1.set_yscale('log')
ax1.legend().set_visible(False)
ax1.set_ylabel('Luminescence, RLU', fontsize=ylabel_size)
ax1.set_xlabel(None)
ax1.tick_params(axis = 'y', labelsize = yticklabel_size)
ax1.set_ylim(4000001, 200000000)
ax1.grid(lw=0.6)

sns.despine(offset=10, trim=False, ax = ax1)

ax1.set_xticks(ticks = np.arange(0, len(sorting_order), 1))
ax1.set_xticklabels(labels = [x for x in sorting_order], fontsize = xticklabel_size)# , rotation=35, ha='right')
# ax1.get_xticklabels()[-1].set_weight('bold')

# statistical annotation
x1, x2 = 0, 1
y, h, col = 8000000, 800000, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data[data["Label"] == "nnLuz WT"]["Sum"].mean(), data[data["Label"] == "nnLuz v4"]["Sum"].mean()) ,1)}-fold\n \
{stat_rounder(mann_whitney_test["nnLuz WT"]["nnLuz v4"])}'

ax1.text((x1+x2)*.5, y-1.5*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

plt.suptitle(r'$\it{Pichia}$' + ' ' + r'$\it{pastoris}$' + f',\n transient expression, ppas_opt\n 50 \u03BCM luciferin treatment (exp1153)', fontsize = title_size)
plt.tight_layout()

plt.savefig(f'', # Insert here path for image to save
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='white')
