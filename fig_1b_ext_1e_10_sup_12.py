# This code with some filereading alterations was used to processing HEK293 data (Figure 1b, Extended Data Figure 1e, 10, Supplementary Figure 12)

import string

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as ss
import scikit_posthocs as scihoc
import seaborn as sns

from os import listdir
from re import search

exp_name = '' # ID of the experiment
path = f'' # Insert here path to raw data file

FBP = {} # Dictionary with plasmid-to-names mapping

def stat_rounder(p):
  if p < 0.0001:
    return 'p = ' + '%.1e' % Decimal(str(p))
  # elif p >= 0.0001 and p < 0.001:
  #   return 'p < 0.001'
  # elif p >= 0.001 and p < 0.01:
  #   return 'p < 0.01'
  # elif p >= 0.01 and p < 0.05:
  #   return 'p < 0.05'
  else:
    return f'p = {str(round(p, 4))}'

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

data = pd.read_csv(f'{path}.csv', encoding= 'unicode_escape') # Insert here path to raw data file
data = data[['Image Number', 'ROI', 'Total Flux [p/s]']]
data['Image Number exact'] = data.apply(lambda row: row['Image Number'].split('_')[-1], axis = 1)

new_data_df = pd.DataFrame(np.nan,
                           index=[data['Image Number exact'].unique()],
                           columns = ['ROI ' + str(x) for x in range(1, len(data['ROI'].unique()) + 1)])

for image_number in data['Image Number exact'].unique():
  for roi in data['ROI'].unique():
    new_data_df.at[image_number, roi] = data[(data['Image Number exact']==image_number) & (data['ROI'] == roi)]['Total Flux [p/s]'].values[0]

new_data_df.columns = [x.replace(' ', '') for x in new_data_df.columns]

roi_names = pd.read_csv(f'{path}.csv', delimiter=';') # Insert here path to roi-to-names mapping
roi_names_dict = {}
for index, row in roi_names.iterrows():
  roi_names_dict[row['ROI']] = row['Set'] + '-' + str(index)

new_data_df.rename(columns = roi_names_dict, inplace = True)
new_data_df.reset_index(inplace = True)
new_data_df.set_index('level_0', inplace = True)

data = new_data_df

time = [x*14 for x in range(len(data.index))] # Depends on exposure time

data = data.transpose()

bg_rows = [x for x in data.index if 'minus' in x]
print(bg_rows)
bg_mean = data[data.index.isin(bg_rows)].mean()
bg_std = data[data.index.isin(bg_rows)].std()
bg = bg_mean - 3 * bg_std
# bg = bg_mean
print(bg)
for index, row in data.iterrows():
  data.loc[index] -= bg

data['Label'] = data.apply(lambda row: row.name.split('-')[0], axis = 1)
data['Sum'] = data.apply(lambda row: row[0:99].sum(), axis = 1)
data.to_csv(f'{path}.csv') # Insert here path for processed data to save


data['Label'] = data.apply(lambda row: FBP[row.Label], axis = 1)
data['Sum_mean'] = data.apply(lambda row: np.divide(row.Sum, len(data_100mm.columns[:-3])), axis = 1) #

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
mann_whitney_test = scihoc.posthoc_conover(data, val_col = 'Sum_mean', group_col = 'Label', p_adjust = 'holm-sidak')


sns.boxplot(data=data,
            x='Label',
            y='Sum_mean', # or Sum for Integral luminescence, RLU
            color = 'white',
            order = sorting_order,
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            ax = ax)

sns.swarmplot(data=data, x='Label', y='Sum_mean', color = 'grey', order = sorting_order, ax = ax) # or Sum for Integral total flux, photons/s

ax.set_yscale('log')
ax.set_xlabel(None)
ax.set_ylabel('Total flux, photons/s', fontsize = ylabel_size)

ax.grid(lw = 0.5)
ax.legend().set_visible(False)

ax.set_xticks(ticks = np.arange(0, len(sorting_order), 1))
ax.set_xticklabels(labels = [x for x in sorting_order])
ax.tick_params(axis = 'y', labelsize = yticklabel_size)
ax.tick_params(axis = 'x', labelsize = xticklabel_size)

y1 = 15000
y2 = 100000
y3 = 300000
y4 = 1000000
y5 = 10000000
ylim_lower = 10000
ylim_higher = 60000000
tick_lower_border = 10

label5 =  r'$\bf{FBP3}$' + '\nmcitHispS\nnnH3H v2\nnnLuz v4\nnnCPH\nNpgA'
label4 = r'$\bf{FBP2}$' + '\nnnHispS\nnnH3H v2\nnnLuz v4\nnnCPH\nNpgA'
label3 = "mcitHispS\nnnH3H v2\nnnLuz WT\nnnCPH\nNpgA"
label2 = "mcitHispS\nnnH3H WT\nnnLuz WT\nnnCPH\nNpgA"
label1 = "nnHispS\nnnH3H WT\nnnLuz WT\nnnCPH\nNpgA"
label0 = r'$\bf{FBP1*}$' + "\nnnHispS\nnnH3H WT\nnnLuz WT\nnnCPH"

ax.set_yticks(ax.get_yticks()[ax.get_yticks() >= tick_lower_border])
ax.set_yticks(ax.get_yticks(minor = True)[ax.get_yticks(minor = True) > tick_lower_border], minor = True)

# statistical annotation
# x1, x2 = 0, 5
# y, h, col = y1, y1/10, 'k'
# ax.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
# p_value = f'{round(np.divide(data[data["Label"] == label5]["Sum"].mean(), data[data["Label"] == label0]["Sum"].mean()) ,1)}-fold, \
# {stat_rounder(mann_whitney_test[label5][label0])}'
# ax.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

x1, x2 = 1, 5
y, h, col = y2, y2/10, 'k'
ax.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data[data["Label"] == label5]["Sum"].mean(), data[data["Label"] == label1]["Sum"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test[label5][label1])}'
ax.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

x1, x2 = 2, 5
y, h, col = y3, y3/10, 'k'
ax.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data[data["Label"] == label5]["Sum"].mean(), data[data["Label"] == label2]["Sum"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test[label5][label2])}'
ax.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

x1, x2 = 3, 5
y, h, col = y4, y4/10, 'k'
ax.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data[data["Label"] == label5]["Sum"].mean(), data[data["Label"] == label3]["Sum"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test[label5][label3])}'
ax.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

x1, x2 = 4, 5
y, h, col = y5, y5/10, 'k'
ax.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data[data["Label"] == label5]["Sum"].mean(), data[data["Label"] == label4]["Sum"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test[label5][label4])}'
ax.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

ax.set_ylim(ylim_lower, ylim_higher)

ns = data.groupby('Label').size()
for label_pos, label in enumerate(sorting_order):
  ax.text(label_pos, ylim_lower, f'N = {ns[label]}', fontsize = xticklabel_size, ha = 'center')

plt.figtext(0.15, 0.8, r'$\bf{Mammalian}$' + '\nHEK293T', fontsize = 24)
sns.despine(offset = 10, trim = False, ax = ax)

plt.savefig(f'{path}.png', # Insert here path for image to save
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='white')

# For kinetics plotting
ax2 = axes['C']

for index, label in enumerate(sorting_order):
  mean = data[data.Label == label].iloc[:, 0:99].mean(axis = 0)
  std = data[data.Label == label].iloc[:, 0:99].std(axis = 0)
  errorfill(x = time, y = mean, yerr = std, color = colors_[index], label = label, ax = ax2)

ax2.legend(loc = 'best', fontsize = legend_size, frameon = False, ncol = 2)
ax2.set_xlabel('Time, min', fontsize = xlabel_size)
ax2.set_ylabel('Total flux, photons/s', fontsize = ylabel_size)


ax2.set_xlim(0)
# ax2.set_ylim(0, 60000000)
ax2.grid(lw = 0.3)
ax2.tick_params(axis = 'both', labelsize = 12)
ticks = np.arange(0, np.max(time), 180)
ax2.set_xticks(ticks = ticks)
ax2.set_xticklabels(labels = [str(int(x/60)) for x in ticks])
