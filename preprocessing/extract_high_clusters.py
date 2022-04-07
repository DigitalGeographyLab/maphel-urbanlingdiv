#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 13:27:19 2021

This script reads data with identified clusters and extracts the high and low clusters.

@author: tuomvais
"""
import glob
import geopandas as gpd
import pandas as pd

# read socioeconomic grid database data in
df = gpd.read_file('RTK_data.gpkg')

# empty list for high and low clusters
highs = []
lows = []

# create empty list for file paths
files = []

# populate list with geopackage file paths to files with bivariate local morans i results
for gpkg in glob.glob('file/path/to/*.gpkg'):
    if 'comb2015' in gpkg:
        files.append(gpkg)
    else:
        pass

# loop over files
for file in files:
    
    # read file
    mdf = gpd.read_file(file)
    
    # get file name
    fn = file.split('/')[-1][:-5]
    
    # extract high shannon/simpson clusters
    mdf = mdf[mdf['sha_sim_cl'] == 1]
    
    # update filename
    fn = fn + '_shasim_high.gpkg'
    
    # add to high list
    if 'comb2015' in fn:
        highs.append(mdf)
    

# loop over files again
for file in files:
    
    # read file
    ldf = gpd.read_file(file)
    
    # get file name
    fn = file.split('/')[-1][:-5]
    
    # extract low clusters
    ldf = ldf[ldf['sha_sim_cl'] == 3]
    
    # update filename
    fn = fn + '_shasim_low.gpkg'
    
    # add to lows list
    lows.append(ldf)

# calculate duplicate geometries    
high_df = pd.concat(highs)
low_df = pd.concat(lows)

# record cluster appearance counts
highcounts = high_df['NRO'].value_counts().rename('stability').reset_index()
lowcounts = low_df['NRO'].value_counts().rename('stability').reset_index()

# drop duplicates
high_df = high_df.drop_duplicates(subset=['NRO'])
low_df = low_df.drop_duplicates(subset=['NRO'])

# join stability series
histab = pd.merge(high_df, highcounts, left_on='NRO', right_on='index')
lowstab = pd.merge(low_df, lowcounts, left_on='NRO', right_on='index')

# save stability to geopackage
histab.to_file('stability_high.gpgk', driver='GPKG')
lowstab.to_file('stability_low.gpgk', driver='GPKG')
