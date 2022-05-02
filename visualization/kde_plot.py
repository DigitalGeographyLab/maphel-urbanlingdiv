#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 13:19:59 2021

This script plots KDE plot of data with cluster information from moran_cluster.py

@author: tuomvais
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Get path to input folder
ap.add_argument("-if", "--inputfolder", required=True,
                help="Path to folder containing files with diversity values with"
                " cluster information.")

# Get path to input folder
ap.add_argument("-o", "--output", required=True,
                help="Path to output graphics file.")

# parse arguments
args = vars(ap.parse_args())

# create empty list for file paths
files = []

# populate list with geopackage file paths
for gpkg in glob.glob(args['inputfolder'] + 'knn8*.gpkg'):
    files.append(gpkg)

# renaming dictionary
cnames = {0:'n.sig.', 1:'High-high', 2:'Low-high', 3:'Low-low', 4:'High-low'}

# color palette dictionary
copal = {'Morning':'goldenrod','Noon':'yellow', 'Afternoon':'darkseagreen', 'Evening':'teal', 'Night':'midnightblue', 'Registry':'gray'}

# Times of day list
titles = ['Afternoon', 'Morning','Night','Noon','Registry','Evening']

# empty list for dataframes
some = []

# empty dataframe for registry data
reg = gpd.GeoDataFrame()

# loop over files
for i, file in enumerate(files):

    # read in
    df = gpd.read_file(file)

    # add column containing times of day
    df['tod'] = titles[i]

    # check where to push df
    if titles[i] == 'Registry':
        # persist registry data
        reg = df
    elif titles[i] != 'Registry':
        # append to social media list
        some.append(df)

# concatenate social media data into one datafrmae
merged = pd.concat(some)

# get sorting order based on times of day
merged.tod = pd.Categorical(merged.tod, categories=['Morning', 'Noon', 'Afternoon', 'Evening', 'Night'])

# sort the data
merged = merged.sort_values('tod')

# rename column to reflect data source
merged = merged.rename(columns={'tod':'Social media'})
reg = reg.rename(columns={'tod':'Registry'})

# reindex to remove duplicate index values
merged = merged.reset_index(drop=True)

# plot kde depicting diversity indices across the datasources and times of day
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

# save the plotted figure as pdf
plt.savefig(args['output'], dpi=300, bbox_inches='tight')
