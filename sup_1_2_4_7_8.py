### This code with was used for processing data for Supplementary Figure 1, 2, 4, 7, 8

import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scipy as sp
import scipy.stats as ss
import scikit_posthocs as scihoc
import seaborn as sns

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
        cbar.ax.tick_params(size=0, labelsize = 15)

        return hax, cbar

path = '' # Insert here path to raw data

replacements = {} # Dictionary with plasmid-to-name mapping

def label_rename(text, replacements = replacements):
  # Create a regular expression from the dictionary keys
  regex = compile("(%s)" % "|".join(map(escape, replacements.keys())))
  # For each match, look-up corresponding value in dictionary
  return regex.sub(lambda mo: replacements[mo.group()], text)


  return label

def stat_rounder(p):
  if p < 0.0001:
    return 'p < 0.0001'
  else:
    return f'p = {str(round(p, 4))}'


data = pd.read_excel('', skiprows = , nrows = ).dropna(axis = 1) # Insert here path for raw data file and coordinates of the needed table
data = data.transpose()

data_to_plot = pd.DataFrame()
for col in data.columns:
  data_to_plot = pd.concat([data_to_plot, data[col].squeeze()])

data_to_plot['cat'] = data_to_plot.apply(lambda row: label_rename(row.name), axis = 1)
data_to_plot.rename(columns={0: 'values'}, inplace=True)
data_to_plot.reset_index(inplace=True)

medianprops = {'color': 'coral'}
capprops = {'color': 'black'}
whiskerprops = {'color': 'black'}
boxprops = {'edgecolor': 'black'}

xlabel_size = 20
ylabel_size = 20
yticklabel_size = 15
xticklabel_size = 15
title_size = 24
suptitle_size = 20
legend_size = 20
signplot_size = 13


fig, axes = plt.subplots(1, 2, figsize = (18, 5), gridspec_kw={'width_ratios': [2.8, 1]})

sorting_order = data_to_plot['cat'].unique()

ax1 = axes[0]
sns.boxplot(data=data_to_plot,
            x='cat',
            y='values',
            color = 'white',
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            order = sorting_order,
            showfliers = False,
            ax = ax1)

sns.swarmplot(data=data_to_plot, x='cat', y='values', order = sorting_order, color = 'grey', ax = ax1)

# ax1.set_yscale('log')
ax1.legend().set_visible(False)
ax1.set_ylabel('Luminescence, RLU', fontsize=ylabel_size)
ax1.set_xlabel(None)
ax1.tick_params(axis = 'y', labelsize = yticklabel_size)

ax1.grid(lw=0.6)

sns.despine(offset=10, trim=False, ax = ax1)

ax1.set_xticks(ticks = np.arange(0, len(sorting_order), 1))
ax1.set_xticklabels(labels = [x for x in sorting_order], fontsize = xticklabel_size, rotation=35, ha='right')
ax1.get_xticklabels()[0].set_weight('bold')

# Example of statistical annotation from Supplementary 2a
mann_whitney_test = scihoc.posthoc_mannwhitney(data_to_plot, val_col = 'values', group_col = 'cat')
# statistical annotation
loc = 30000
x1, x2 = 0, 1
y, h, col = loc, loc/10, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data_to_plot[data_to_plot["cat"] == "nnLuz v3"]["values"].mean(), data_to_plot[data_to_plot["cat"] == "nnLuz WT"]["values"].mean()) ,1)}-fold\n \
{stat_rounder(mann_whitney_test["nnLuz v3"]["nnLuz WT"])}'
ax1.text((x1+x2)*.5, y-2*h, p_value, ha='center', va='top', fontsize = 18)


ax2 = axes[1]

############ <><><><><><><><>< STATS
q1 = data_to_plot.copy()
q1['datapoint']=q1['values'].apply(float)

data_ = [q1.loc[ids, 'datapoint'].values for ids in q1.groupby('cat').groups.values()]
H, p = ss.kruskal(*data_)
print(f"Results of Kruskal, p={p}, H={H}")
res = scihoc.posthoc_conover(q1, val_col='datapoint', group_col='cat', p_adjust = 'holm-sidak')
sorting_order2 = sorting_order

res=res.reindex(sorting_order2).transpose().reindex(sorting_order2)
heatmap_args = {'linewidths': 0.25, 'linecolor': '0.5', 'clip_on': False, 'square': True, 'cbar_ax_bbox': [0.92, 0.20, 0.02, 0.3]}
ax2 = sign_plot(res,**heatmap_args, ax = ax2, ticksize = signplot_size)



plt.subplots_adjust(hspace=1.2, wspace = 0.4, top = 0.8)
plt.suptitle(r'$\it{Escherichia}$' + ' ' + r'$\it{coli}$' + ', lysates,\nafter 100 \u03BCM luciferin treatment', fontsize = title_size) # Insert here a title. An example is provided

plt.savefig('', # Insert here path for figure to save
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='white')
