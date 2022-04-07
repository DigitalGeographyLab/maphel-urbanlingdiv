#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 13:19:59 2021

@author: tuomvais
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob

# create empty list for file paths
files = []

# populate list with geopackage file paths
for gpkg in glob.glob('/home/tuomvais/GIS/maphel_thirdplace/combined_some/autocorr/knn8*.gpkg'):
    files.append(gpkg)

# renaming dictionary
cnames = {0:'n.sig.', 1:'High-high', 2:'Low-high', 3:'Low-low', 4:'High-low'}

# color palette dictionary
copal = {'Morning':'goldenrod','Noon':'yellow', 'Afternoon':'darkseagreen', 'Evening':'teal', 'Night':'midnightblue', 'Registry':'gray'}

# Times of day list
titles = ['Afternoon', 'Morning','Night','Noon','Registry','Evening']

# empty list for dataframes
some = []
reg = gpd.GeoDataFrame()

# loop over files 
for i, file in enumerate(files):
    
    # read in
    df = gpd.read_file(file)
    
    # add column
    df['tod'] = titles[i]
    
    # check where to push df
    if titles[i] == 'Registry':
        reg = df
    elif titles[i] != 'Registry':
        # append to social media list 
        some.append(df)

# merge together social media data
merged = pd.concat(some)

# get sorting order
merged.tod = pd.Categorical(merged.tod, categories=['Morning', 'Noon', 'Afternoon', 'Evening', 'Night'])

# sort
merged = merged.sort_values('tod')

# rename
merged = merged.rename(columns={'tod':'Social media'})
reg = reg.rename(columns={'tod':'Registry'})

# reindex to remove duplicate index values
merged = merged.reset_index(drop=True)

# plot kde
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

sns.set()
fig, axes = plt.subplots(1, 2, figsize=(15,6), sharex=True, sharey=True)
fig.suptitle('Language diversity across times of day')
reg_line = mlines.Line2D([], [], color='grey', linestyle='--', label='Registry data')
sns.kdeplot(ax=axes[0], data=merged, x='shannon_scaled', hue='Social media',
            common_norm=False, palette='Spectral', alpha=.8, lw=2.2,)
sns.kdeplot(ax=axes[0], data=reg, x='shannon_scaled', color='black',
            common_norm=False, linestyle='--',alpha=.5, lw=2.2,)
plt.legend(loc='upper right')
plt.legend(loc='right', bbox_to_anchor=(-0.2, 0.65), handles=[reg_line])
axes[0].set_title('a. Scaled Shannon entropy')
axes[0].set_xlabel('')
sns.kdeplot(ax=axes[1], data=merged, x='simpson', hue='Social media',
            common_norm=False, palette='Spectral', alpha=.8, lw=2.2, legend=False)
sns.kdeplot(ax=axes[1], data=reg, x='simpson', color='black',
            common_norm=False, linestyle='--',alpha=.5, lw=2.2,)
axes[1].set_title('b. Simpson diversity')
axes[1].set_xlabel('')
#axes[1].get_legend().remove()
plt.savefig('/home/tuomvais/GIS/maphel_thirdplace/plots/kde_diversities_timeofday_2.pdf', dpi=300, bbox_inches='tight')

ax = sns.kdeplot(data=merged, x='shannon_scaled', hue='Time of day', palette='rocket_r', alpha=.6, lw=2, bw_adjust=.6, multiple='layer')

sns.kdeplot(data=merged, x='shannon', hue='tod', fill=True, alpha=.5, common_norm=False, bw=.1, multiple='fill', )




sns.set()
fig, axes = plt.subplots(1, 2, figsize=(15,6), sharex=True, sharey=True)
fig.suptitle('Language diversity across times of day')
reg_patch = mpatches.Patch(color='black', alpha=0.3, label='Registry data')
sns.histplot(ax=axes[0], data=merged, x='shannon_scaled', hue='Social media',
            palette='Spectral', alpha=.8, stat='probability', common_norm=False,
            multiple='stack', element='step', bins=33)
sns.histplot(ax=axes[0], data=reg, x='shannon_scaled', color='black',
            alpha=.3, stat='probability', element='step', bins=33)
plt.legend(loc='upper right')
plt.legend(loc='right', bbox_to_anchor=(-0.2, 0.65), handles=[reg_patch])
axes[0].set_title('a. Scaled Shannon entropy')
axes[0].set_xlabel('')
sns.histplot(ax=axes[1], data=merged, x='simpson', hue='Social media',
            palette='Spectral', alpha=.8, stat='probability', common_norm=False,
            multiple='stack', element='step', bins=33, legend=False)
sns.histplot(ax=axes[1], data=reg, x='simpson', color='black',
            alpha=.3, stat='probability', element='step', bins=33)
axes[1].set_title('b. Simpson diversity')
axes[1].set_xlabel('')
#axes[1].get_legend().remove()
plt.savefig('/home/tuomvais/GIS/maphel_thirdplace/plots/hist_diversities_timeofday.pdf', dpi=300, bbox_inches='tight')



sns.set()
fig, axes = plt.subplots(1, 2, figsize=(15,6), sharex=True, sharey=True)
fig.suptitle('Language diversity across times of day')
reg_patch = mpatches.Patch(color='black', alpha=0.3, label='Registry data')
sns.histplot(ax=axes[0], data=merged, x='shannon_scaled', hue='Social media',
            palette='Spectral', alpha=.8, stat='probability', common_norm=False,
            multiple='fill', element='step', bins=33)
sns.histplot(ax=axes[0], data=reg, x='shannon_scaled', color='black',
            alpha=.3, stat='probability', element='step', bins=33)
plt.legend(loc='upper right')
plt.legend(loc='right', bbox_to_anchor=(-0.2, 0.65), handles=[reg_patch])
axes[0].set_title('a. Scaled Shannon entropy')
axes[0].set_xlabel('')
sns.histplot(ax=axes[1], data=merged, x='simpson', hue='Social media',
            palette='Spectral', alpha=.8, stat='probability', common_norm=False,
            multiple='fill', element='step', bins=33, legend=False)
sns.histplot(ax=axes[1], data=reg, x='simpson', color='black',
            alpha=.3, stat='probability', element='step', bins=33)
axes[1].set_title('b. Simpson diversity')
axes[1].set_xlabel('')
#axes[1].get_legend().remove()
plt.savefig('/home/tuomvais/GIS/maphel_thirdplace/plots/hist_fill_diversities_timeofday.pdf', dpi=300, bbox_inches='tight')