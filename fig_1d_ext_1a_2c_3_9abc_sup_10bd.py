# This code with some filereading alterations was used for processing and plotting data from BY-2, obtainded using TECAN (Figure 1d, Extended Data 1a, 2c, 9abc, Supplementary 10bd)

import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
import scipy.stats as ss
import scikit_posthocs as scihoc
from decimal import Decimal
from re import search

from typing import Union, List, Tuple

from matplotlib import colors
from matplotlib.axes import SubplotBase
from matplotlib.colorbar import ColorbarBase, Colorbar
from matplotlib.colors import ListedColormap

def sign_plot(
        x: Union[List, np.ndarray, pd.DataFrame],
        g: Union[List, np.ndarray] = None,
        flat: bool = False,
        labels: bool = True,
        cmap: List = None,
        cbar_ax_bbox: List = None,
        ax: SubplotBase = None,
        ticksize = 12,
        **kwargs) -> Union[SubplotBase, Tuple[SubplotBase, Colorbar]]:
    """Significance plot, a heatmap of p values (based on Seaborn).
    Parameters
    ----------
    x : Union[List, np.ndarray, DataFrame]
        If flat is False (default), x must be an array, any object exposing
        the array interface, containing p values. If flat is True, x must be
        a sign_array (returned by :py:meth:`scikit_posthocs.sign_array`
        function).
    g : Union[List, np.ndarray]
        An array, any object exposing the array interface, containing
        group names.
    flat : bool
        If `flat` is True, plots a significance array as a heatmap using
        seaborn. If `flat` is False (default), plots an array of p values.
        Non-flat mode is useful if you need to  differentiate significance
        levels visually. It is the preferred mode.
    labels : bool
        Plot axes labels (default) or not.
    cmap : list
        1) If flat is False (default):
        List consisting of five elements, that will be exported to
        ListedColormap method of matplotlib. First is for diagonal
        elements, second is for non-significant elements, third is for
        p < 0.001, fourth is for p < 0.01, fifth is for p < 0.05.
        2) If flat is True:
        List consisting of three elements, that will be exported to
        ListedColormap method of matplotlib. First is for diagonal
        elements, second is for non-significant elements, third is for
        significant ones.
        3) If not defined, default colormaps will be used.
    cbar_ax_bbox : list
        Colorbar axes position rect [left, bottom, width, height] where
        all quantities are in fractions of figure width and height.
        Refer to `matplotlib.figure.Figure.add_axes` for more information.
        Default is [0.95, 0.35, 0.04, 0.3].
    ax : SubplotBase
        Axes in which to draw the plot, otherwise use the currently-active
        Axes.
    kwargs
        Keyword arguments to be passed to seaborn heatmap method. These
        keyword args cannot be used: cbar, vmin, vmax, center.
    Returns
    -------
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes object with the heatmap.
    cbar : matplotlib.colorbar.Colorbar
        ColorBar object if `flat` is set to False.
    Examples
    --------
    >>> x = np.array([[ 1, 1, 1],
                      [ 1, 1, 0],
                      [ 1, 0, 1]])
    >>> ph.sign_plot(x, flat = True)
    """
    for key in ['cbar', 'vmin', 'vmax', 'center']:
        if key in kwargs:
            del kwargs[key]

    if isinstance(x, pd.DataFrame):
        df = x.copy()
    else:
        x = np.array(x)
        g = g or np.arange(x.shape[0])
        df = pd.DataFrame(np.copy(x), index=g, columns=g)

    dtype = df.values.dtype

    if not np.issubdtype(dtype, np.integer) and flat:
        raise ValueError("X should be a sign_array or DataFrame of integers")
    elif not np.issubdtype(dtype, np.floating) and not flat:
        raise ValueError("X should be an array or DataFrame of float p values")

    if not cmap and flat:
        # format: diagonal, non-significant, significant
        cmap = ['1', '#fbd7d4', '#1a9641']
    elif not cmap and not flat:
        # format: diagonal, non-significant, p<0.001, p<0.01, p<0.05
        cmap = ['1', '#fbd7d4', '#005a32', '#238b45', '#a1d99b']

    if flat:
        np.fill_diagonal(df.values, -1)
        hax = sns.heatmap(df, vmin=-1, vmax=1, cmap=ListedColormap(cmap),
                      cbar=False, ax=ax, **kwargs)
        if not labels:
            hax.set_xlabel('')
            hax.set_ylabel('')
        return hax

    else:
        df[(x < 0.001) & (x >= 0)] = 1
        df[(x < 0.01) & (x >= 0.001)] = 2
        df[(x < 0.05) & (x >= 0.01)] = 3
        df[(x >= 0.05)] = 0

        np.fill_diagonal(df.values, -1)

        if len(cmap) != 5:
            raise ValueError("Cmap list must contain 5 items")

        hax = sns.heatmap(
            df, vmin=-1, vmax=3, cmap=ListedColormap(cmap), center=1,
            cbar=False, ax=ax, **kwargs)
        if not labels:
            hax.set_xlabel('')
            hax.set_ylabel('')
        hax.tick_params(axis = 'both', labelsize = ticksize)

        cbar_ax = hax.figure.add_axes(cbar_ax_bbox or [0.95, 0.35, 0.04, 0.3])
        cbar = ColorbarBase(cbar_ax, cmap=(ListedColormap(cmap[2:] + [cmap[1]])), norm=colors.NoNorm(),
                            boundaries=[0, 1, 2, 3, 4])
        cbar.set_ticks(list(np.linspace(0, 3, 4)))
        cbar.set_ticklabels(['p < 0.001', 'p < 0.01', 'p < 0.05', 'NS'])
        cbar.outline.set_linewidth(1)
        cbar.outline.set_edgecolor('0.5')
        cbar.ax.tick_params(size=0)

        return hax, cbar

path = '' # Insert here path to raw data
exp = '' # Insert here experiment ID

FBP = {} # Dictionary with plasmid-to-name mapping

def stat_rounder(p):
  if p < 0.0001:
    return 'p < 0.0001'
  else:
    return f'p={str(round(p, 4))}'

def detect_table(data):
    start = None
    for index, val in enumerate(data.iloc[:, 0]):
        if val == '<>' and data.iloc[index + 1, 1] != 'None':
            start = index

    if start:
        data_r = data.iloc[start + 1:start + 9, :13]
        data_r.columns = data.iloc[start]
        data_r = data_r.set_index('<>')

        return data_r


def plate_process(overall_table, filename, plate):
  overall_table['values'] = np.nan

  df = pd.read_excel(f'{path}{filename}')
  df = detect_table(df)
  for row_index, row in overall_table[overall_table.Plate == plate].iterrows():
      row_name = row.Row
      col_name = str(row.Col)
      overall_table.iloc[row_index, -1] = df[df.index == row_name][col_name]
  try:
    assert len(overall_table[overall_table['values'] == np.nan].index.values) == 0
  except AssertionError:
    nans = overall_table[overall_table['values'] == np.nan]
    nans_names = nans.Label.unique()
    print(f'In {"|".join(nans_names)} some values were not provided.')
  to_plot = overall_table
  to_plot = to_plot[to_plot['Ignore'] == False]
  return to_plot

def get_gly(text):
  if '(' in text:
    text = text.replace(')', '')
    gly = text.split('(')[0]
  else:
    gly = ''
  return gly

def get_label(text):
  text = text.split('_')[0]
  if '(' in text:
    text = text.replace(')', '')
    label = text.split('(')[1]
  else:
    label = text
  return label

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

overall_table = pd.read_excel('') # Insert here path to table with info about plasmids in each plate well, generated from private Airtable Base
data_to_plot = pd.DataFrame()
parts = [] # Listed all part of data, if more than one
for part in parts:
  data_to_plot1 = plate_process(overall_table, '', plate = part) # Insert here path to raw data file
  data_to_plot1 = data_to_plot1[data_to_plot1.Plate == part]
  data_to_plot = pd.concat([data_to_plot, data_to_plot2])

data_to_plot['Label for plotting (from AC)'] = data_to_plot.apply(lambda row: row['Label for plotting (from AC)'].replace('OD', ''), axis = 1)
data_to_plot['Label'] = data_to_plot.apply(lambda row: get_label(row['Label for plotting (from AC)']), axis = 1)
data_to_plot['Label'] = data_to_plot.apply(lambda row: FBP[row['Label']], axis = 1)
data_to_plot['Dilution'] = data_to_plot.apply(lambda row: row['Label for plotting (from AC)'].split('_')[-1], axis = 1)
data_to_plot['Gly'] = data_to_plot.apply(lambda row: get_gly(row['Label for plotting (from AC)']), axis = 1)
data_to_plot.to_csv('') # Insert here path for processed data to save


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

label3 = r'$\bf{FBP3}$' + '\nmcitHispS\nNpgA\nnnLuz v4\nnnCPH\nnnH3H v2'
label2 = r'$\bf{FBP2}$' + '\nnnHispS\nNpgA\nnnLuz v4\nnnCPH\nnnH3H v2'
label1 = '\nnnHispS\n\nnnLuz v4\nnnCPH\nnnH3H v2'
label0 = r'$\bf{FBP1}$' + '\nnnHispS\n\nnnLuz WT\nnnCPH\nnnH3H WT'

dil = # selected OD of agrocombinations, if there are several options
data_ = data_to_plot[(data_to_plot.Dilution == dil)]
data_['Label'] = data_.apply(lambda row: get_label(row['Label for plotting (from AC)']), axis = 1)
data_['Label'] = data_.apply(lambda row: FBP[row['Label']], axis = 1)


mann_whitney_test = scihoc.posthoc_mannwhitney(data_, val_col = 'values', group_col = 'Label', p_adjust = 'holm-sidak')
sorting_order_all = data_.groupby('Label').median().sort_values('values', ascending = True).index.values

fig, axes = plt.subplots(1, figsize = (12, 8))
ax1 = axes
sns.boxplot(data=data_,
            x='Label',
            y='values',
            color = 'white',
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            order = sorting_order_all,
            showfliers = False,
            ax = ax1)

sns.swarmplot(data=data_, x='Label', y='values', color = 'grey', order = sorting_order_all, ax = ax1)

ax1.set_yscale('log')
ax1.legend().set_visible(False)
ax1.set_ylabel('Luminescence, RLU', fontsize=ylabel_size)
ax1.set_xlabel(None)
ax1.tick_params(axis = 'y', labelsize = yticklabel_size)

ax1.grid(lw=0.6)

sns.despine(offset=10, trim=False, ax = ax1)

ax1.set_xticks(ticks = np.arange(0, len(sorting_order_all), 1))
ax1.set_xticklabels(labels = [x for x in sorting_order_all], fontsize = xticklabel_size)

y1 = 60000
y2 = 200000
y3 = 1000000
ylim_lower = 22000
ylim_higher = 80000000
tick_lower_border = 10



ax1.set_yticks(ax1.get_yticks()[ax1.get_yticks() >= tick_lower_border])
ax1.set_yticks(ax1.get_yticks(minor = True)[ax1.get_yticks(minor = True) > tick_lower_border], minor = True)

# Example of statistical annotation from Figure 1d
# statistical annotation
x1, x2 = 0, 3
y, h, col = y1, y1/10, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data_[data_["Label"] == label3]["values"].mean(), data_[data_["Label"] == label0]["values"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test[label3][label0])}'
ax1.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

x1, x2 = 1, 3
y, h, col = y2, y2/10, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data_[data_["Label"] == label3]["values"].mean(), data_[data_["Label"] == label1]["values"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test[label3][label1])}'
ax1.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

x1, x2 = 2, 3
y, h, col = y3, y3/10, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data_[data_["Label"] == label3]["values"].mean(), data_[data_["Label"] == label2]["values"].mean()) ,1)}-fold, \
{stat_rounder(mann_whitney_test[label3][label2])}'
ax1.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

ns = data.groupby('Label').size()
for label_pos, label in enumerate(sorting_order_all):
  ax1.text(label_pos, ylim_lower, f'N = {ns[label]}', fontsize = xticklabel_size, ha = 'center')

ax1.set_ylim(ylim_lower, ylim_higher)


plt.figtext(0.15, 0.8, r'$\bf{Plant}$' + '\nBY-2', fontsize = 24) # Insert here a title. Example is provided

plt.savefig('', # Insert here path for figure to save
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='white')




# In case of plotting sign-plot with Conover statistics
 ax2 = axes[1]

############ <><><><><><><><>< STATS
q1 = data_.copy()
q1['datapoint']=q1['values'].apply(float)

data_ = [q1.loc[ids, 'datapoint'].values for ids in q1.groupby('Label').groups.values()]
H, p = ss.kruskal(*data_)
print(f"Results of Kruskal, p={p}, H={H}")
res = scihoc.posthoc_conover(q1, val_col='datapoint', group_col='Label', p_adjust = 'holm-sidak')
sorting_order2 = sorting_order

res=res.reindex(sorting_order2).transpose().reindex(sorting_order2)
heatmap_args = {'linewidths': 0.25, 'linecolor': '0.5', 'clip_on': False, 'square': True, 'cbar_ax_bbox': [0.92, 0.20, 0.021, 0.28]}
ax2 = sign_plot(res,**heatmap_args, ax = ax2, ticksize = signplot_size)
