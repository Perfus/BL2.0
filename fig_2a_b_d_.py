### This code with some filereading alterations was used for processing data captured on Sony Alpha ILCE-7M3 camera (Figure 2a, 2b, 2d, )

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

selected = ['pX037', 'pNK497']
data_selected = data[(data.group.isin(selected)) & (data.iso == 'iso400')]
data_selected.to_csv('/content/gdrive/Shareddrives/✨ NOUKA ✨/Ноучные эксперименты/nouka26. Paper about bioluminescent system 2.0/01. Preparing draft/Data/data15/exp948/exp948_data_selected.csv')
