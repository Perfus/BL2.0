### This code was used to remake plant images in Extended Data 6a, Supplementaty Figures 14a, 19ac, 20a

import numpy as np
import matplotlib.pyplot as plt
import tifffile as tf
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

vmin = 0
linthresh = 50 # Depends on background values of photo

ax1 = axes[0]
ax1pos = ax1.get_position()
ax1left = ax1pos.x0
ax1bottom = ax1pos.y0
ax1width = ax1pos.width
ax1height = ax1pos.height

photo1_ratio = ax1width / ax1height
ax1.set_position([0, ax1bottom - 0.15, ax1width + ax1left, ax1height + ax1left / photo1_ratio])

fn = '' # Insert path to tif file here
fname = f'{path}{fn}'
im = tf.imread(fname)
im = im.astype(float)
im -=2048
im = im.clip(0)
vmax = np.percentile(im.ravel(),99.999) # Or set it manually

image = ax1.imshow(im, cmap='inferno', norm=colors.SymLogNorm(linthresh=linthresh, vmin=vmin, vmax=vmax, base = 10), aspect='equal')
ax1.axis('off')
ax1.plot([230, 900], [850, 850], color = 'white', linewidth = 3)
ax1.plot([1250, 2000], [850, 850], color = 'white', linewidth = 3)
ax1.text(445, 950, 'FBP2', color = 'white', fontsize = 20)
ax1.text(1525, 950, 'FBP3', color = 'white', fontsize = 20)

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="10%", pad = 0.5)
cbar = plt.colorbar(image, orientation='horizontal', extend='both', shrink = 0.3, cax = cax)
cax.tick_params(labelsize=14)
cbar.set_label('Luminescence, RLU', fontsize = 20)
