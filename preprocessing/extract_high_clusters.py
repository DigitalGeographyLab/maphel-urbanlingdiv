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
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Get path to input file
ap.add_argument("-i", "--input", required=True,
                help="Path to input RTK database file (type geopackage).")

# Get path to input file
ap.add_argument("-if", "--inputfolder", required=True,
                help="Path to folder with geopackagee files containing cluster"
                " information. Example: /path/to/folder/ ")

# Get path to output file
ap.add_argument("-of", "--outputfolder", required=True,
                help="Path to output folder Example: /path/to/outputfolder/.")

# parse arguments
args = vars(ap.parse_args())

# read socioeconomic grid database data in
df = gpd.read_file(args['input'])

# empty list for high and low clusters
highs = []
lows = []

# create empty list for file paths
files = []

# populate list with geopackage file paths to files with bivariate local morans i results
for gpkg in glob.glob(args['inputfolder'] + '*.gpkg'):
    files.append(gpkg)

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
histab.to_file(args['outputfolder'] + 'stability_high.gpgk', driver='GPKG')
lowstab.to_file(args['outputfolder'] + 'stability_low.gpgk', driver='GPKG')
