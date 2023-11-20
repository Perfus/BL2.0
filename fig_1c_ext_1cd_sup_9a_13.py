### This code with some filereading alterations was used for processing data from full plants captured on Sony Alpha ILCE-7M3 camera (Figure 1c, Extended Data Figure 1cd, Supplementary Figure 9a, 13)

import string

from decimal import Decimal
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as ss
import scikit_posthocs as scihoc
import seaborn as sns

from os import listdir
from re import search

exp_name = '' # ID of the experiment
path = f'' # Insert here path to the raw data files

def errorfill(x, y, yerr, color=None, alpha_fill=0.3, ax=None, label=None):
    ax = ax if ax is not None else plt.gca()
    if color is None:
        color = ax._get_lines.color_cycle.next()
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin, ymax = yerr
    ax.plot(x, y, color=color, label=label)
    ax.fill_between(x, ymax, ymin, color=color, alpha=alpha_fill)

colors_ = ['darkblue',
          'darkgreen',
          'magenta',
          'teal',
          'tan',
          'mediumaquamarine',
          'lightcoral',
          'gray',
          'gold',
          'deepskyblue']

FBP = {} # Dictionary with plasmid-to-name mapping

def stat_rounder(p):
  return 'p = ' + '%.1e' % Decimal(str(p))
  # elif p >= 0.0001 and p < 0.001:
  #   return 'p < 0.001'
  # elif p >= 0.001 and p < 0.01:
  #   return 'p < 0.01'
  # elif p >= 0.01 and p < 0.05:
  # #   return 'p < 0.05'
  # else:
  #   return f'p = {str(round(p, 4))}'

annotations = pd.read_csv(f'{path}.csv', delimiter = ';', index_col = 'Order') # Insert here path to file with roi-to-names
anno_dict = dict((str(index), f"{row['Name'].replace('  ', ' ')}-{row['ppas']}-{row['Spot']}") for index, row in annotations.iterrows())


data = pd.read_csv(f'{path}.csv', index_col = ' ') # Insert here path to raw data file
data = data.transpose()
data = data.loc[:, :31]
time = [x*3 for x in range(len(data.columns))]
data['long_label'] = data.apply(lambda row: row.name.replace('Mean', ''), axis = 1)
data['long_label'] = data.apply(lambda row: anno_dict[row['long_label']], axis = 1)
data.set_index('long_label', inplace = True)
bg_rows = [x for x in data.index if 'BG' in x]

bg_mean = data[data.index.isin(bg_rows)].mean()
bg_std = data[data.index.isin(bg_rows)].std()
bg = bg_mean
# bg = bg_mean - 3 * bg_std # This BG was used when some samples were not luminescent and some datapoints was <= mean BG roi values, i.e. Extended Data 1cd
print(bg)
for index, row in data.iterrows():
  data_100mm.loc[index] -= bg

data = data[~data.index.isin(bg_rows)]
data['Label'] = data.apply(lambda row: row.name.split('-')[0], axis = 1)
data['ppas'] = data.apply(lambda row: row.name.split('-')[1], axis = 1)
data['Sum'] = data.apply(lambda row: row[0:31].sum(), axis = 1) # Insert selected number of frames instead of 31 depending on kinetics
data['Sum_norm'] = data.apply(lambda row: row['Sum']/len(row[0:31]), axis = 1)
data['Max'] = data.apply(lambda row: row[0:31].max(), axis = 1)

data['Label'] = data.apply(lambda row: FBP[row.Label], axis = 1)

medianprops = {'color': 'coral'}
capprops = {'color': 'black'}
whiskerprops = {'color': 'black'}
boxprops = {'edgecolor': 'black'}

xlabel_size = 20
ylabel_size = 22
yticklabel_size = 15
xticklabel_size = 18
title_size = 24
suptitle_size = 20
legend_size = 20
signplot_size = 16
letter_size = 24

fig, axes = plt.subplots(1, figsize = (12, 8))
ax = axes

sorting_order = data.groupby('Label').median().sort_values('Sum', ascending = True).index.values
mann_whitney_test = scihoc.posthoc_mannwhitney(data, val_col = 'Sum_norm', group_col = 'Label', p_adjust = 'holm-sidak')

sns.boxplot(data=data,
            x='Label',
            y='Sum_norm', # Or Sum for Integral luminescence, RLU
            color = 'white',
            order = sorting_order,
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            ax = ax)

sns.swarmplot(data=data, x='Label', y='Sum_norm', hue = 'ppas', order = sorting_order, ax = ax) # Or Sum for Integral luminescence, RLU

ax.set_yscale('log')
ax.set_xlabel(None)
ax.set_ylabel('Luminescence, RLU', fontsize = ylabel_size)


ax.grid(lw = 0.5)
ax.legend().set_visible(False)

ax.set_xticks(ticks = np.arange(0, len(sorting_order), 1))
ax.set_xticklabels(labels = [x for x in sorting_order])
ax.tick_params(axis = 'y', labelsize = yticklabel_size)
ax.tick_params(axis = 'x', labelsize = xticklabel_size)

y1 = 2
y2 = 20
y3 = 200
y4 = 1000
ylim_lower = 1
ylim_higher = 30000
tick_lower_border = 0

label4 = r"$\bf{FBP3}$" + "\nmcitHispS\nnnH3H v2\nnnLuz v4\nNpgA\n"
label3 = "mcitHispS\nnnH3H v2\nnnLuz WT\nNpgA"
label2 = "mcitHispS\nnnH3H WT\nnnLuz WT\nNpgA"
label1 = "nnHispS\nnnH3H WT\nnnLuz WT\nNpgA"
label0 = r"$\bf{FBP1*}$" + "\nnnHispS\nnnH3H WT\nnnLuz WT"

ax.set_yticks(ax.get_yticks()[ax.get_yticks() >= tick_lower_border])
ax.set_yticks(ax.get_yticks(minor = True)[ax.get_yticks(minor = True) > tick_lower_border], minor = True)

# Example of statistical annotations from Figure 1c
# statistical annotation
# x1, x2 = 0, 4
# y, h, col = y1, y1/10, 'k'
# ax.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
# p_value = f'{round(np.divide(data[data["Label"] == label4]["Sum"].mean(), data[data["Label"] == label0]["Sum"].mean()) ,1)}-fold, \
# {stat_rounder(mann_whitney_test[label4][label0])}'
# ax.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

x1, x2 = 1, 4
y, h, col = y2, y2/10, 'k'
ax.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data[data["Label"] == label4]["Max"].mean(), data[data["Label"] == label1]["Max"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test[label4][label1])}'
ax.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

x1, x2 = 2, 4
y, h, col = y3, y3/10, 'k'
ax.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data[data["Label"] == label4]["Max"].mean(), data[data["Label"] == label2]["Max"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test[label4][label2])}'
ax.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

x1, x2 = 3, 4
y, h, col = y4, y4/10, 'k'
ax.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data[data["Label"] == label4]["Max"].mean(), data[data["Label"] == label3]["Max"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test[label4][label3])}'
ax.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

ax.set_ylim(ylim_lower, ylim_higher)

ns = data_100mm.groupby('Label').size()
for label_pos, label in enumerate(sorting_order):
  ax.text(label_pos, ylim_lower, f'N = {ns[label]}', fontsize = xticklabel_size, ha = 'center')

plt.figtext(0.15, 0.8, r'$\bf{Yeast}$' + '\n' + r'$\it{Pichia}$' + ' ' + r'$\it{pastoris}$', fontsize = 24) # Insert here a title. Example is provided
sns.despine(offset = 10, trim = False, ax = ax)


plt.savefig(f'', # Insert here path for image to save
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='white')


# If kinetics were ploted

ax2 = axes['C']

for index, label in enumerate(sorting_order):
  mean = data[data.Label == label].iloc[:, 0:31].mean(axis = 0)
  std = data[data.Label == label].iloc[:, 0:31].std(axis = 0)
  errorfill(x = time, y = mean, yerr = std, color = colors_[index], label = label, ax = ax2)

ax2.legend(loc = 'best', fontsize = legend_size, frameon = False, bbox_to_anchor=[1.01, 1.], mode = 'expand')
ax2.set_xlabel('Time, min', fontsize = xlabel_size)
ax2.set_ylabel('Luminescence, RLU', fontsize = ylabel_size)


ax2.set_xlim(0)
ax2.grid(lw = 0.3)
ax2.tick_params(axis = 'both', labelsize = 12)
# ticks = np.arange(0, np.max(time), 180)
# ax2.set_xticks(ticks = ticks)
# ax2.set_xticklabels(labels = [str(int(x/60)) for x in ticks])
