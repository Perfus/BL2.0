### This code was used for processing data from full plants captured on Sony Alpha ILCE-7M3 camera (Supplementary Figure 14b, 19cd)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss

from os import listdir
from re import search

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

FBP = {} # Dictionary with plasmid-to-name mapping

filelist = listdir(f'{path}/')
filelist = [x for x in filelist if '.csv' in x]
full = pd.DataFrame()
for file in filelist:
  data = pd.read_csv(f'{path}/{file}', usecols = ['Label', 'Area', 'RawIntDen'])
  data['Line'] = data.apply(lambda row: row.Label.split('-')[0].split(':')[1], axis = 1)

  bgs = []
  for i, row in data[data.Line == 'bg'].iterrows():
    bgs.append(row['RawIntDen']/row['Area'])
  bg = np.mean(bgs)

  data['bg_sub'] = data.apply(lambda row: row.RawIntDen - (row.Area * bg), axis = 1)
  data = data[data.Line != 'bg']
  data['plasmid'] = data.apply(lambda row: FBP[group(row.Line.upper())], axis = 1)
  data['group'] = data.apply(lambda row: group(row.Line.upper()), axis = 1)
  full = pd.concat([full, data])
full.to_csv('') # Insert here path for processed data to save

xlabel_size = 20
ylabel_size = 20
yticklabel_size = 15
xticklabel_size = 18
title_size = 24
suptitle_size = 20
legend_size = 20
signplot_size = 16
letter_size = 24

fig, axes = plt.subplots(1, figsize = (10, 6))
medianprops = {'color': 'coral'}
capprops = {'color': 'black'}
whiskerprops = {'color': 'black'}
boxprops = {'edgecolor': 'black'}

ax1 = axes
sns.boxplot(data=full,
            x='group',
            y='bg_sub',
            color = 'white',
            order = selected,
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            showfliers = False,
            ax = ax1)

sns.swarmplot(data=full, x='group', y='bg_sub', hue = 'Line', order = selected, ax = ax1)

ax1.set_yscale('log')
ax1.set_xlabel(None)
ax1.set_ylabel(u'Luminescense, RLU', fontsize = 20)
ax1.legend().set_visible(False)

ax1.grid(lw = 0.3)


ax1.set_xticks(ticks = np.arange(0, len(selected), 1))
ax1.set_xticklabels(labels = [FBP[x] for x in selected])
ax1.tick_params(axis = 'y', labelsize = 14)
ax1.tick_params(axis = 'x', labelsize = 18)
ax1.set_title(r'$\it{Nicotiana}$' + ' ' + r'$\it{benthamiana}$' + ', transgenic lines', ha='center', va='center', fontsize = 22) # Insert here a title. Example is provided
ax1.set_ylim(100001)

# Example of fold annotation from Supplementary Figure 1b
# statistical annotation
x1, x2 = 0, 1
y, h, col = 400000, 40000, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(full[full["group"] == "pNK511"]["bg_sub"].mean(), full[full["group"] == "pX037"]["bg_sub"].mean()) ,1)}-fold' \
#({stat_rounder(mann_whitney_test_y["pNK497"]["pX037"])})'
ax1.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

# statistical annotation
x1, x2 = 0, 2
y, h, col = 250000, 25000, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(full[full["group"] == "pNK497"]["bg_sub"].mean(), full[full["group"] == "pX037"]["bg_sub"].mean()) ,1)}-fold' \
#({stat_rounder(mann_whitney_test_y["pNK497"]["pX037"])})'
ax1.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

# statistical annotation
x1, x2 = 0, 3
y, h, col = 150000, 15000, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(full[full["group"] == "pNK3074"]["bg_sub"].mean(), full[full["group"] == "pX037"]["bg_sub"].mean()) ,1)}-fold' \
#({stat_rounder(mann_whitney_test_y["pNK3074"]["pX037"])})'
ax1.text((x1+x2)*.5, y-1.8*h, p_value, ha='center', va='top', fontsize = xticklabel_size)

sns.despine(offset = 10, trim = False, ax = ax1)

plt.savefig('', # Insert here path for figure to save
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='white')


 
