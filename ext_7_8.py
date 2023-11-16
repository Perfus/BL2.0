### This code with some filereading alterations was used for processing data of FBP3 vs NanoLuc/FFLuc comparison using Tecan (Extended Data 7, 8)

import string
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D
from scipy import integrate
import scipy.stats as ss
import scikit_posthocs as scihoc
import seaborn as sns

from os import listdir
from re import search

exp_name = '' # Insert here ID of experiment
path = f'/{exp_name}/{exp_name}/' # Insert here path to raw data

def stat_rounder(p):
  if p < 0.0001:
    return 'p < 0.0001'
  else:
    return f'p={str(round(p, 4))}'

def errorfill(x, y, yerr, color=None, alpha_fill=0.3, ax=None, label=None):
    ax = ax if ax is not None else plt.gca()
    if color is None:
        color = ax._get_lines.color_cycle.next()
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin, ymax = yerr
    ax.plot(x, y, color=color, label=label, linewidth = 2)
    ax.fill_between(x, ymax, ymin, color=color, alpha=alpha_fill)

colors_ = ['darkblue',
          'darkgreen',
          'magenta',
          'teal',
          'tan',
          'mediumaquamarine',
          'lightcoral',
          'gray']

# Creating dictionary with wells-to-names mapping for discrete file with raw data or for all raw data
order = pd.read_excel('') # Insert here path to file with wells-to-names mapping
dict_names = {}
for well in order['Well']:
  dict_names[well] = order[order.Well == well]['Description'].values[0]

def labeler(name):
  if 'nanoLuc' in name:
    return 'NanoLuc'
  elif 'Time' in name:
    return 'Time'
  elif 'FBP3_MS_medium' in name:
    return r'$\bf{FBP3}$' + '\nmcitHispS\nnnH3H v2\nnnLuz v4\nnnCPH\nNpgA'
  else:
    return 'Ignore'


# Example of raw data from 1 file processing
nrows = # Insert here the quantity of time points in raw data file
plate_x = pd.read_excel('', skiprows = 49, nrows = nrows, index_col = 'Cycle Nr.') # Insert here path to file with raw data
time_x = plate_x['Time [s]']
plate_x_dict = {}
plate_x_dict['Time'] = plate_x_dict
for col in plate_x.columns[2:]:
  if not math.isnan(plate_x[col][1]):
    plate_nluc_MS_dict[col] = plate_x[col]
plate_x_dict = pd.DataFrame(plate_x_dict)
plate_x_dict = plate_x_dict.rename(columns = dict_names)
plate_x_dict = plate_x_dict.transpose()
plate_x_dict['Label'] = plate_x_dict.apply(lambda row: labeler(row.name), axis = 1)
plate_x_dict = plate_x_dict[plate_nluc_MS.Label != 'Ignore']

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

plates = [] # List all plate objects obtained using previous block of code
target_time = 1800 # Time for signal integrating
data_box = pd.DataFrame()

for plate in plates:
  time_point = find_nearest(plate.values[0].tolist()[:-1], target_time))
  print(f'{time_point} in {plate}'
  plate['NPSum20min'] = plate.apply(lambda row: np.trapz(y=row[:time_point], x = plate_other_100.values[0].tolist()[:time_point]), axis = 1)
  for_sns_plate = plate[['NPSum20min', 'Label']]
  for_sns_plate.set_axis(['NPSum20min', 'Label'], axis=1, inplace=True)
  data_box = pd.concat([data_box, for_sns_plate])

data_box = data_box[['minus' not in s for s in data_box.index]]

mann_whitney_test = scihoc.posthoc_mannwhitney(data_box, val_col = 'NPSum20min', group_col = 'Label', p_adjust = 'holm-sidak')

fig, axes = plt.subplots(1, 2, figsize = (18, 8))

sorting_order = [] # Manually set list of sample names

medianprops = {'color': 'coral'}
capprops = {'color': 'black'}
whiskerprops = {'color': 'black'}
boxprops = {'edgecolor': 'black'}

ax1 = axes[0]
sns.boxplot(data=data_box,
            x='Label',
            y='NPSum20min',
            color = 'white',
            order = sorting_order,
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            ax = ax1)

sns.swarmplot(data=data_box, x='Label', y='NPSum20min', color='grey', order = sorting_order, ax = ax1)

ax1.set_yscale('log')
ax1.set_xlabel(None)
ax1.set_ylabel('Integral luminescence, RLU', fontsize = 20)
ax1.grid(lw = 0.3)

ax1.set_xticks(ticks = np.arange(0, len(sorting_order), 1))
ax1.set_xticklabels(labels = [x for x in sorting_order])
ax1.tick_params(axis = 'y', labelsize = 15)
ax1.tick_params(axis = 'x', labelsize = 18)

FBP3 = r'$\bf{FBP3}$' + '\nmcitHispS\nnnH3H v2\nnnLuz v4\nnnCPH\nNpgA'
# Example of statistical annotation in Extended Data Figure 7
x1, x2 = 0, 2
y, h, col = 40000000000, 6000000000, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data_box[data_box["Label"] == "NanoLuc"]["NPSum20min"].mean(), data_box[data_box["Label"] == FBP3]["NPSum20min"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test["NanoLuc"][FBP3])}'
ax1.text((x1+x2)*.5, y-1.5*h, p_value, ha='center', va='top', fontsize = 18)

x1, x2 = 1, 2
y, h, col = 1500000000, 200000000, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data_box[data_box["Label"] == FBP3]["NPSum20min"].mean(), data_box[data_box["Label"] == "FFLuc"]["NPSum20min"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test["FFLuc"][FBP3])}'
ax1.text((x1+x2)*.5, y-1.5*h, p_value, ha='center', va='top', fontsize = 18)


ax2 = axes[1]

# Example of plotting one errorfill plot
nanoluc_timepoint = # Insert here a corresponding time point for required time 
for index, label in enumerate(plate_x['Label'].unique()):
  if 'NanoLuc' in label:
    mean = plate_x[plate_nluc_MS.Label == label].iloc[:, 0:nanoluc_timepoint].mean(axis = 0)
    std = plate_x[plate_nluc_MS.Label == label].iloc[:, 0:nanoluc_timepoint].std(axis = 0)
    label_ = label
    color_ = 'darkblue'
    errorfill(x = plate_x.iloc[0, :nanoluc_timepoint].values, y = mean, yerr = std, color = color_, label = label_, ax = ax2)

ax2.legend(loc = 1, fontsize = 18, frameon = False, ncol = 2)
ax2.set_xlabel('Time, min', fontsize = 20)
ax2.set_ylabel('Luminescence, RLU', fontsize = 20)

ax2.set_yscale('log', base = 10)

ax2.set_xlim(0)
ax2.grid(lw = 0.3)
ax2.tick_params(axis = 'both', labelsize = 15)
ticks = ax2.get_xticks()
ax2.set_xticklabels(labels = [str(int(x/60)) for x in ticks])

axes = axes.flat
for n, ax in enumerate(axes):
    ax.text(-0.2, 1.1,
            string.ascii_lowercase[n],
            transform=ax.transAxes,
            size=35)

plt.subplots_adjust(hspace=0.5)
plt.suptitle('', fontsize = 22) # Insert here title of the plot

sns.despine(offset = 10, trim = False)
plt.tight_layout(pad = 2.0, h_pad = 0.5)

plt.savefig('', # Insert here path for figure to save 
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='white')


