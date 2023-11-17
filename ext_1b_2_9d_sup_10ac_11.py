# This code was used for processing and plotting data from plant leaves captured by Sony Alpha ILCE-7M3 camera (Extended Data 1b, 2, 9d, Supplementary 10ac, 11)

import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss
import scikit_posthocs as scihoc

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
background_columns = ('BG1', 'BG2', 'BG3')

FBP = {} # Dictionary with plasmid-to-name mapping

def stat_rounder(p):
  if p < 0.0001:
    return 'p < 0.0001'
  else:
    return f'p={str(round(p, 4))}'

def subtract_background(leaf_analysis_table,
                        leaf_analysis_ID,
                        data):

    columns = [c for c in data.columns if 'Mean' in c]

    leaf_ID = leaf_analysis_table['Layout'][leaf_analysis_ID]
    leaves_per_photo = int(leaf_analysis_table['Leaves/photo'][leaf_analysis_ID])
    leaves = [f'{leaf_ID}-{leaf_number+1}' for leaf_number in range(0, leaves_per_photo)]
    print('Leaves in the file:', leaves)

    roi_order = leaf_analysis_table['Spot ROI order'][leaf_analysis_ID]
    spot_columns = [f'Label {int(spot)}' for spot in roi_order.split('-')]
    print('Spot order:', spot_columns)

    all_conditions = leaf_analysis_table.loc[leaf_analysis_ID][spot_columns].values
    conditions_exp = [c for c in all_conditions if c != 'nothing']
    print('Conditions: ', conditions_exp)

    # Background ROI order is assumed to be standard: (L * Spots) - (L * 3* BG_leaf) - 3* BG_out
    conditions_bg = ['BG1', 'BG2', 'BG3']
    columns_bg_out = columns[-3:]

    # assumes ROI were selected in this specific order in ImageJ
    conditions_names = []
    for leaf in leaves:
        for cond in conditions_exp:
            conditions_names.append(f'{leaf} {cond}')
    for leaf in leaves:
        for cond in conditions_bg:
            conditions_names.append(f'{leaf} {cond}')
    print(len(columns))
    print(len(conditions_names))
    assert len(columns) == len(conditions_names) + 3


    # renaming columns
    column_names_dict = dict(zip(columns, conditions_names))
    data.rename(columns=column_names_dict, inplace=True)

    # subtracting mean value of background outside of leaves
    for c in data[column_names_dict.values()]:
        bg = data[columns_bg_out].mean(axis=1) - 3 * data[columns_bg_out].std(axis=1)
        if data[c][1] < bg[1]:
          print(c)
          print(data[c][1])
          print(bg[1])
        data[c] = data[c] - (data[columns_bg_out].mean(axis=1) - 3 * data[columns_bg_out].std(axis=1))

    # removing bg_out and other columns we won't need
    data = data[column_names_dict.values()]

    ID_leaf_replicate_condition = [(leaf_analysis_ID,
                                    c.split()[0].split('-')[0],
                                    c.split()[0].split('-')[1],
                                    ' '.join(c.split()[1:])) for c in data.columns]

    multiindex = pd.MultiIndex.from_tuples(ID_leaf_replicate_condition, names=['analysis_ID',
                                                                               'leaf_layout',
                                                                               'replicate',
                                                                               'infiltration'])
    data.columns = multiindex


    return data


def df_transpose(data):
  data = data.transpose()
  data.columns = [f'Mean{c}' for c in list(data.columns)]
  data = pd.DataFrame(data.loc['Mean']).transpose()
  return data

def data_prep(data, overall_table, name, ac_labels = None, norm = None):
  leaf_layouts_dropped = []
  for analysis_id in overall_table.index.values:

      leaf_layout_id = overall_table.loc[analysis_id]['Layout']

      if leaf_layout_id not in leaf_layouts_dropped:

          leaf_replicates = set([replicate for _, _, replicate, _ in data.loc(axis=1)[analysis_id, leaf_layout_id, :, :].columns])
          print(leaf_replicates)
          infiltration_columns = [infiltration for _, _, _, infiltration in data.loc(axis=1)[analysis_id, leaf_layout_id, :, :].columns]

  # dropping background columns
  for c in background_columns:
      data.drop(c, axis=1, level=3, inplace=True)
  data.columns = data.columns.remove_unused_levels()

  to_save = data.transpose()
  to_save = to_save.reset_index()
  if ac_labels:
    to_save['label'] = to_save.apply(lambda row: ac_labels[row['infiltration']], axis = 1)
  to_save.to_hdf(f'{path}{name}_not_normed.hdf5', 'data')
  return to_save

overall_table = pd.read_excel('', index_col = 'ID', usecols=lambda x: 'Unnamed' not in x) # Insert here path to table with info about spots, generated from private Airtable Base
data = pd.DataFrame()

for leaf_analysis_ID in overall_table.index:
    url = f"{path}{leaf_analysis_ID}_{overall_table['ImageJ analysis csv'][leaf_analysis_ID]}.csv"
    raw_data = pd.read_csv(url, index_col = ' ')
    new_data = subtract_background(overall_table,
                                   leaf_analysis_ID,
                                   raw_data)
    data = pd.concat([data, new_data], axis=1)
data.to_hdf(f'{path}{exp}data__all_leaves.hdf5', key='data', mode='w')

pre_data = pd.read_hdf(f'{path}{exp}data__all_leaves.hdf5', 'data')
data = data_prep(data = pre_data, overall_table = overall_table, name = exp)
data['OD'] = data.apply(lambda row: row['infiltration'].split('_')[-1].replace('OD', ''), axis = 1)
data['label'] = data.apply(lambda row: FBP[row['infiltration'].split('_')[0]], axis = 1)


xlabel_size = 20
ylabel_size = 20
yticklabel_size = 15
xticklabel_size = 18
fold_size = 18
title_size = 22
suptitle_size = 24
legend_size = 20
signplot_size = 16
letter_size = 40

medianprops = {'color': 'coral'}
capprops = {'color': 'black'}
whiskerprops = {'color': 'black'}
boxprops = {'edgecolor': 'black'}

selected = [] # List of plasmids to be considered on the plot
dil = # selected OD of agrocombinations, if there are several options

data_ = data[(data.OD == dil) & (data.label.isin(selected))]
fig, axes = plt.subplots(1, 2, figsize = (19, 8))

sorting_order = data_.groupby('label').median().sort_values('Mean', ascending = True).index.values

medianprops = {'color': 'coral'}
capprops = {'color': 'black'}
whiskerprops = {'color': 'black'}
boxprops = {'edgecolor': 'black'}

ax1 = axes[0]
sns.boxplot(data=data_,
            x='label',
            y='Mean',
            color = 'white',
            order = sorting_order,
            medianprops = medianprops,
            capprops = capprops,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            showfliers = False,
            ax = ax1)

sns.swarmplot(data=data_, x='label', y='Mean', color='grey', order = sorting_order, ax = ax1)

ax1.set_yscale('log')
ax1.set_xlabel(None)
ax1.set_ylabel('Luminescence, RLU', fontsize = 20)

ax1.grid(lw = 0.3)


ax1.set_xticks(ticks = np.arange(0, len(sorting_order), 1))
ax1.set_xticklabels(labels = [x for x in sorting_order])
ax1.set_yticks(ax1.get_yticks()[2:])
ax1.set_yticks(ax1.get_yticks(minor = True)[ax1.get_yticks(minor = True) > 1], minor = True)
ax1.tick_params(axis = 'y', labelsize = 15)
ax1.tick_params(axis = 'x', labelsize = 18)
ax1.set_ylim(0.0400000001, 2000)


# Example of statistical annotation from Supplementary Figure 11a
pNK497 = r'$\bf{FBP2}$' + '\nnnHispS\nnnH3H v2\nnnLuz v4\nnnCPH\nNpgA'
pNK3071 = r'$\bf{FBP3}$' + '\nmcitHispS\nnnH3H v2\nnnLuz v4\nnnCPH\nNpgA'
mann_whitney_test = scihoc.posthoc_mannwhitney(data_, val_col = Mean, group_col = 'label', p_adjust = 'holm-sidak')
# statistical annotation
x1, x2 = 0, 1
y, h, col = 40, 4, 'k'
ax1.plot([x1, x1, x2, x2], [y, y-h, y-h, y], lw=1.5, c=col)
p_value = f'{round(np.divide(data_[data_["label"] == pNK3071][1].median(), data_[data_["label"] == pNK497][1].median()) ,1)}-fold\n \
{stat_rounder(mann_whitney_test[pNK3071][pNK497])}'

ax1.text((x1+x2)*.5, y-1.5*h, p_value, ha='center', va='top', fontsize = 13)


ax2 = axes[1]
# Example of sign plot from Conover test from Supplementary Figure 10a
########### <><><><><><><><>< STATS
q1 = data_.copy()
q1['datapoint']=q1['Mean'].apply(float)
q1['name'] = q1.apply(lambda row: row['label'], axis = 1)

data = [q1.loc[ids, 'datapoint'].values for ids in q1.groupby('name').groups.values()]
H, p = ss.kruskal(*data)
print(f"Results of Kruskal, p={p}, H={H}")
res = scihoc.posthoc_conover(q1, val_col='datapoint', group_col='name', p_adjust = 'holm-sidak')
sorting_order2 = sorting_order
sorting_order2
res=res.reindex(sorting_order2).transpose().reindex(sorting_order2)
heatmap_args = {'linewidths': 0.25, 'linecolor': '0.5', 'clip_on': False, 'square': True, 'cbar_ax_bbox': [0.85, 0.69, 0.027, 0.13]}
hax, cbar = sign_plot(res,**heatmap_args, ax = ax2, ticksize = 16)
hax.tick_params('y', labelrotation = 0)


  plt.savefig('', # Insert here path for figure to save
              dpi = 400,
              bbox_inches='tight',
              transparent=False,
              facecolor='white')

