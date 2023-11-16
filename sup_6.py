### This code was used to fit data to Michaelis-Menten equation and determine Vmax and Km parameters for nnLuz WT and nnLuz v4

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns

from scipy.optimize import minimize, curve_fit

exp_name = '' # Insert here ID of experiment
path_raw = f'/{exp_name}/{exp_name} raw data/' # Insert here path to raw data

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

def plate_process(filename,  anno_dict):
  print(filename)
  df = pd.read_excel(f'{path_raw}{filename}')
  df = detect_table(df)
  df = df.dropna(axis=1, how='all')
  df = df.loc[:, df.columns != 'Temp. [°C]']

  df = df.transpose()

  df['Label'] = df.apply(lambda row: ''.join(anno_dict[anno_dict.Order == row.name]['Name'].values).replace('_ppas_opt', '').replace('_', ' ').replace('wt', 'WT'), axis = 1)
  df['ppas'] = df.apply(lambda row: ''.join(anno_dict[anno_dict.Order == row.name]['ppas'].values), axis = 1)
  df['use'] = df.apply(lambda row: anno_dict[anno_dict.Order == row.name]['use'].values, axis = 1)
  df.loc[df.Label=='none','Label']="ingore 0 ignore"

  return df

# MM model
def v(s, v_max, k_m):
  return (v_max * s) / (k_m + s)

anno_dict = pd.read_excel('') # Insert here path to file with wells-to-names mapping

filelist = [x for x in os.listdir(path_raw) if 'part' in x]
ppas_data = pd.DataFrame()
for file in filelist:
  data = plate_process(file,
                       anno_dict)
  ppas_data = pd.concat([ppas_data, data[data.columns[:57].tolist() + ['Label', 'ppas','use']].iloc[1:, :]], axis = 0) # 57 frames is the first 15 minutes of data collecting
time = pd.DataFrame(data.iloc[0, :-3]).transpose()

ppas_data['Sum'] = ppas_data.apply(lambda row: row[:-3].sum(), axis = 1)
ppas_data['conc'] = ppas_data.apply(lambda row: float(row.Label.split(' ')[-2]), axis = 1)
ppas_data['use'] = ppas_data['use'].apply(lambda x: x[0])
ppas_data.sort_values(by = ['ppas', 'conc'], inplace = True)
ppas_data = pd.concat([time, ppas_data])
ppas_data.to_csv(f'{path_raw}{exp_name}_data.csv')


anno_dict_hibit = pd.read_excel('') # Insert here path to file with wells-to-names mapping of HiBit data file

data_hibit = plate_process('', # Insert here path to file HiBit values
                          anno_dict_hibit)

data_hibit = data_hibit[data_hibit.columns[:41].tolist() + ['Label', 'ppas']] # 41 frames is the first20 minutes of HiBit data collecting
bg_mean = data_hibit[data_hibit.Label == 'minus'].mean()
bg = bg_mean
for index, row in data_hibit.iloc[1:, :].iterrows():
  data_hibit.loc[index] -= bg

data_hibit['Label'] = data_hibit.apply(lambda row: ''.join(anno_dict_hibit[anno_dict_hibit.Order == row.name]['Name'].values).replace('_ppas_opt', '').replace('_', ' ').replace('wt', 'WT'), axis = 1)
data_hibit['ppas'] = data_hibit.apply(lambda row: ''.join(anno_dict_hibit[anno_dict_hibit.Order == row.name]['ppas'].values), axis = 1)
data_hibit = data_hibit[data_hibit.Label != 'minus']
data_hibit['Sum'] = data_hibit.apply(lambda row: row[:-2].sum(), axis = 1)
data_hibit.to_csv(f'{path_raw}{exp_name}_data_hibit.csv')

to_km = ppas_data.query('use>0 and conc <50').iloc[0:, :].copy() # Choose concentration below 50 mkM
for index, row in to_km.iterrows():
  to_km.loc[index, to_km.columns[:-5]] /= hibit[hibit.index == row.ppas]['Sum'].values[0]
to_km['Sum'] = to_km.apply(lambda row: row[:-5].sum(), axis = 1)
to_km['Luz'] = to_km.apply(lambda row: ' '.join(row.Label.split(' ')[0:2]), axis = 1)
to_plot = to_km[['Sum', 'conc', 'Luz']].groupby(['Luz', 'conc']).describe()

xlabel_size = 22
ylabel_size = 22
yticklabel_size = 18
xticklabel_size = 12
title_size = 24
suptitle_size = 22
legend_size = 20
signplot_size = 13


plt.figure(figsize=(7,7))
ax = plt.subplot()
for index, label in enumerate(to_plot.index.get_level_values(0).unique().values):
  mean = to_plot[to_plot.index.get_level_values(0) == label]['Sum', 'mean'].values
  std = to_plot[to_plot.index.get_level_values(0) == label]['Sum', 'std'].values
  conc = to_plot.index.get_level_values(1).unique().values

  s_real = to_km[to_km.Luz == label][['conc']].values
  v_real = to_km[to_km.Luz == label][['Sum']].values
  popt, pcov = curve_fit(v, s_real.ravel(), v_real.ravel(), p0=[1, 1])

  # Calculate the standard deviations (sqrt of the diagonals of the covariance matrix)
  perr = np.sqrt(np.diag(pcov))
  # Print the results
  print(f"curve_fit() v_max = {popt[0]:.3f} ± {perr[0]:.3f}")
  print(f"curve_fit() k_m = {popt[1]:.3f} ± {perr[1]:.3f}")

  s_plot = np.linspace(0, 30, 100)
  ax.plot(s_plot, v(s_plot, res.x[0], res.x[1]), color = colors_[index])

  plot = ax.errorbar(x = conc, y = mean, yerr = std, color = colors_[index], label = label + f':\nVmax = {popt[0]:.5f} ± {perr[0]:.5f}\
  \nKm = {popt[1]:.2f} ± {perr[1]:.2f} \u03BCM\n',
                      lw = 2,fmt='.', marker = 's', markersize = 5 , capsize = 9)
  plot[-1][0].set_linestyle('--')

ax.legend(loc = 2, fontsize = legend_size, frameon = False,bbox_to_anchor=[1,1.05])
ax.set_xlabel('Concentration, \u03BCM', fontsize = xlabel_size)
ax.set_ylabel('Normalized\n integral luminescence, RLU', fontsize = ylabel_size)

ax.set_xscale('linear')
ax.set_xlim(0)
ax.set_ylim(0)
ax.grid(lw = 0.5)
ax.tick_params(axis = 'both', labelsize = 12)
sns.despine(offset=10, trim=False, ax = ax)
plt.suptitle(r'$\it{Pichia}$' + ' ' + r'$\it{pastoris}$' + f', lysates', fontsize = title_size)

plt.savefig('', # Insert here path for figure to save
            dpi = 400,
            bbox_inches='tight',
            transparent=False,
            facecolor='white')
