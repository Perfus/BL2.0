# This code was used for processing and plotting data from IVIS (Figure 2c, Supplementary 15, 16, 17)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss
import scikit_posthocs as scihoc
import string

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

xlabel_size = 20
ylabel_size = 20
yticklabel_size = 15
xticklabel_size = 17
title_size = 24
suptitle_size = 20
legend_size = 20
signplot_size = 16
letter_size = 24

path = '' # Insert here path to raw data

to_group = {'pX037': ['NB021'],
            'Background': ['BG', 'bg'],
           } # Dictionary with roi-to-plasmid mapping

fbp = {'pX037' : r'$\bf{FBP1}$',
       'pNK3074': r'$\bf{FBP3}$',
       'pNK511': 'nnHispS\nnnH3H v2\nnnLuz v4\nnnCPH',
       'pNK497': r'$\bf{FBP2}$'} # Dictionary with plasmid-to-name mapping

def group(row):
  for k, v in to_group.items():
    if row in v:
      return k

selected = ['pX037', 'pNK511', 'pNK497', 'pNK3074'] # List of plasmids to be considered on the plot

chosel = ['part1', 'part2', 'part3', 'part4', 'part5']
data = pd.DataFrame()
for index_, file in enumerate(chosel):
  new_data = pd.read_csv(f'{path}/_{file}.csv', delimiter = ';', decimal = ',', usecols=['Avg Radiance [p/s/cm?/sr]', 'Line']) # Insert here common name of group

  # Background subtracting
  bg  = new_data[new_data.Line == 'bg']['Avg Radiance [p/s/cm?/sr]'].mean()

  for index, row in new_data.iterrows():
    new_data['Avg Radiance [p/s/cm?/sr]'].loc[index] -= bg

  new_data = new_data[new_data.Line != 'bg']

  norm_rows = [x for x in new_data.index if new_data['Line'].loc[x] == 'nb021_L']
  norm_value = new_data['Avg Radiance [p/s/cm?/sr]'].loc[norm_rows].mean()

  # If normalizing
  for row in new_data.index:
    new_data['Avg Radiance [p/s/cm?/sr]'].loc[row] = np.divide(new_data['Avg Radiance [p/s/cm?/sr]'].loc[row], norm_value)
  new_data['Line + group'] = new_data.apply(lambda row: row['Line'].split('_')[0].upper() + '_' + file, axis = 1)

  if index_ > 0:
    new_data = new_data[~new_data.Line.isin(['nb021_L'])]

  data = pd.concat([data, new_data], axis = 0, ignore_index = True)
data['FBP'] = data.apply(lambda row: group(row['Line'].split('_')[0].upper()), axis = 1)
data = data[data.group.isin(selected)]

mann_whitney_test = scihoc.posthoc_mannwhitney(leaves, val_col = 'Avg Radiance [p/s/cm?/sr]', group_col = 'group', p_adjust = 'holm-sidak')


fig, axes = plt.subplots(1, 2 figsize = (18, 8))

sorting_order = data.groupby('group').median().sort_values('Avg Radiance [p/s/cm?/sr]', ascending = True).index.values
ax3 = axes[0]
sns.boxplot(data=data,
            x='FBP',
            y='Avg Radiance [p/s/cm?/sr]',
            color = 'white',
            order = sorting_order,
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            showfliers = False,
            ax = ax3)

sns.swarmplot(data=data, x='group', y='Avg Radiance [p/s/cm?/sr]', hue = 'Line', order = sorting_order, ax = ax3)

ax3.set_yscale('log')
ax3.set_xlabel(None)
ax3.set_ylabel(u'Normalized average radiance, p/s/cm\u00B2/sr', fontsize = ylabel_size)
ax3.legend().set_visible(False)

ax3.grid(lw = 0.3)


ax3.set_xticks(ticks = np.arange(0, len(sorting_order), 1))
ax3.set_xticklabels(labels = [FBP[x] for x in sorting_order])
ax3.tick_params(axis = 'y', labelsize = yticklabel_size)
ax3.tick_params(axis = 'x', labelsize = xticklabel_size)

# Example of statistical annotation from Supplementary Figure 17
# statistical annotation
x1, x2 = 0, 1
y, h, col = 3000000, 400000, 'k'
ax3.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = str(round(np.divide(data[data["FBP"] == r"$\bf{FBP2}$" + "\nnnHispS\nnnH3H v2\nnnLuz v4\nnnCPH\nNpgA"]["Avg Radiance [p/s/cm?/sr]"].mean(), data[data["FBP"] == r"$\bf{FBP1}$" + "\nnnHispS\nnnH3H WT\nnnLuz WT\nnnCPH"]["Avg Radiance [p/s/cm?/sr]"].mean()) ,1))+'-fold\n \
'+stat_rounder(mann_whitney_test[r"$\bf{FBP2}$" + "\nnnHispS\nnnH3H v2\nnnLuz v4\nnnCPH\nNpgA"][r"$\bf{FBP1}$" + "\nnnHispS\nnnH3H WT\nnnLuz WT\nnnCPH"])
ax3.text((x1+x2)*.5, y-5*h, p_value, ha='center', va='bottom', fontsize = 18)
ax3.set_ylim(800000)

ax4 = axes[1]
########### <><><><><><><><>< STATS

q1 = data.copy()
q1['datapoint']=q1['Avg Radiance [p/s/cm?/sr]'].apply(float)
q1['name'] = q1.apply(lambda row: fbp[row['FBP']], axis = 1)

data_stat = [q1.loc[ids, 'datapoint'].values for ids in q1.groupby('group').groups.values()]
H, p = ss.kruskal(*data_stat)
print(f"Results of Kruskal, p={p}, H={H}")
res = scihoc.posthoc_conover(q1, val_col='datapoint', group_col='name', p_adjust = 'holm-sidak')
sorting_order2 = [fbp[x] for x in sorting_order]
sorting_order2
res=res.reindex(sorting_order2).transpose().reindex(sorting_order2)
heatmap_args = {'linewidths': 0.25, 'linecolor': '0.5', 'clip_on': False, 'square': True, 'cbar_ax_bbox': [0.999, 0.35, 0.045, 0.3]}
hax2, cbar2 = sign_plot(res,**heatmap_args, ax = ax4, ticksize = signplot_size)
hax2.tick_params('y', labelrotation = 0)
# Comment previous block and re-plot if Kruskal-Wallis H0 hypothesis is not rejected


plt.subplots_adjust(hspace=3)
plt.suptitle(r'$\it{Nicotiana}$' + ' ' + r'$\it{benthamiana}$' + ', transgenic lines', fontsize = 24) # Insert here a title. Example is provided
sns.despine(offset = 10, trim = False, ax = ax3)
plt.tight_layout(pad = 3)
plt.savefig('', dpi = 900, bbox_inches='tight', transparent=False, # Insert here path for figure to save
            facecolor='white')
