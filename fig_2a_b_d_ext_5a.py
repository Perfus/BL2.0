### This code with some filereading alterations was used for processing data captured on Sony Alpha ILCE-7M3 camera (Figure 2a, 2b, 2d, Extended Data Figure 5a, )

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scikit_posthocs as scihoc

path = '' # Insert here path to raw data

to_group = {} # Dictionary with roi-to-plasmid mapping

FBP = {} # Dictionary with plasmid-to-name mapping

def group(row):
  for k, v in to_group.items():
    if row in v:
      return k

filelist = [x for x in listdir(path) if '.csv' in x and '.png' not in x]
data = pd.DataFrame()

for file in filelist:

  new_data = pd.read_csv(f'{path}{file}', index_col = ' ', usecols = [' ', 'Label', 'Mean'])
  new_data['nLabel'] = new_data.apply(lambda row: row['Label'].split(':')[1].split('-')[0], axis = 1)
  new_data['iso'] = file.split('_')[1]

  data_rows = [x for x in new_data.index if 'bg' not in new_data['nLabel'].loc[x]]
  bg_rows = [x for x in new_data.index if 'bg' in new_data['nLabel'].loc[x]]
  bg = new_data[new_data.index.isin(bg_rows)].mean()

  for i in range(len(new_data.index)):
    new_data['Mean'].iloc[i] -= bg

  new_data = new_data[new_data.index.isin(data_rows)]

  data = pd.concat([data, new_data], ignore_index = True)

data['group'] = data.apply(lambda row: group(row['nLabel']), axis = 1)

selected = [] # List of plasmids to plot
data_selected = data[(data.group.isin(selected)) & (data.iso == 'iso400')]
data_selected.to_csv('') # Insert here path to save processed data

mann_whitney_test = scihoc.posthoc_mannwhitney(data_selected, val_col = 'Mean', group_col = 'group', p_adjust = 'holm-sidak')

fig, axes = plt.subplots(1, figsize = (18, 5))

medianprops = {'color': 'coral'}
capprops = {'color': 'white'}
whiskerprops = {'color': 'white'}
boxprops = {'edgecolor': 'white'}

sorting_order = data_selected.groupby('group').median().sort_values('Mean', ascending=True).index.values

ax2 = axes[1]

sns.boxplot(data=data_selected,
            x='group',
            y='Mean',
            color = 'black',
            order = sorting_order,
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            ax = ax2)
sns.swarmplot(data=data_selected, x='group', y='Mean', hue = 'nLabel', order = sorting_order, ax = ax2)

ax2.set_yscale('log')
ax2.set_xlabel(None)
ax2.set_ylabel(u'Luminescence, RLU', fontsize = ylabel_size)
ax2.set_title(r'$\it{Nicotiana}$' + ' ' + r'$\it{tabacum}$' + '\n', fontsize = title_size, color = '#ffdfff')  # Insert here a title. Example is provided

ax2.legend().set_visible(False)

ax2.grid(lw = 0.3)
ax2.set_ylim(0.7)
ax2.set_xticks(ticks = np.arange(0, len(sorting_order), 1))
ax2.set_xticklabels(labels = [FBP[x] for x in sorting_order])
ax2.tick_params(axis = 'y', labelsize = yticklabel_size)
ax2.tick_params(axis = 'x', labelsize = xticklabel_size)

# Example of statistical annotation from Figure 2b
# statistical annotation
x1, x2 = 0, 1
y, h, col = 1.3, 0.2, 'white'
ax2.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data_selected[data_selected["group"] == "pNK497"]["Mean"].mean(), data_selected[data_selected["group"] == "pX037"]["Mean"].mean()) ,1)}-fold\n \
{stat_rounder(mann_whitney_test["pNK497"]["pX037"])}'
ax2.text((x1+x2)*.5, y-2.1*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

plt.subplots_adjust(hspace=4)
sns.despine(offset = 10, trim = False)
plt.tight_layout(w_pad = 2.9)

plt.savefig('', # Insert here path for figure to save
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='black')
