#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 12:10:25 2021

This script calculates grid cell diversities for the whole social media data without separating times of day.

@author: waeiski
"""
import geopandas as gpd
import pandas as pd
import pysal
import libpysal
import skbio.diversity.alpha as sk
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Get path to register input file
ap.add_argument("-r", "--register", required=True,
                help="Path to input register file (type geopackage).")

# Get path to grid database input file
ap.add_argument("-i", "--input", required=True,
                help="Path to combined Twitter and Instagram point feature file"
                " (type geopackage).")

# Get path to output file
ap.add_argument("-o", "--output", required=True,
                help="Path to output file (type geopackage).")

# parse arguments
args = vars(ap.parse_args())

# function to scale values
def scale_minmax(series):
    '''
    Parameters
    ----------
    series : pandas series containing lang counts in a list.
    Returns
    -------
    Min-Max scaled values
    '''
    # get values
    series = series.values
    # reshape for scikit
    series = series.reshape(-1,1)
    # scale
    scaled = MinMaxScaler().fit_transform(series)
    return scaled

# read grid in
grid = gpd.read_file(args['register'])

# reduce column clutter
grid = grid[['NRO', 'KUNTA', 'KUNTANRO', 'geometry']]

# convert grid id to integer
grid['NRO'] = grid['NRO'].astype(int)

# get contiguity
gW = libpysal.weights.Queen.from_dataframe(grid, idVariable='NRO')

# read timeofday df in
df = gpd.read_file(args['input'])

# join with grid
df = gpd.sjoin(df, grid, op='intersects')

# get grid ids
idlist = df['NRO'].value_counts().index.tolist()

# set language dictionary
langdict = {}

# loop over grid ids
for gid in idlist:
    # get neighboring grid ids of current grid
    neigh_ids = list(gW[gid].keys())
    # extend list with current grid id
    neigh_ids.extend([gid])
    # fetch neighborhood grid cells for current grid id into dataframe
    ndf = df[df['NRO'].isin(neigh_ids)]
    # extract languages used in neighborhood grid cells
    langs = list(zip(ndf['language'].value_counts().index,
                     ndf['language'].value_counts().values))
    # add language use to language dictionary
    langdict[gid] = langs

# generate a language use dataframe
langdf = pd.DataFrame(langdict.items(), columns=['grid_id', 'langs'])

# loop over language dataframe
print('[INFO] - Calculating diversity indices..')
for i, row in langdf.iterrows():
    # get language counts from tuples
    counts = [x[1] for x in row['langs']]
    # calculate diversities
    langdf.at[i, 'sents'] = sum(counts)
    langdf.at[i, 'unique'] = sk.observed_otus(counts)
    langdf.at[i, 'singletons'] = sk.singles(counts)
    langdf.at[i, 'berger'] = sk.berger_parker_d(counts)
    langdf.at[i, 'dominance'] = sk.dominance(counts)
    langdf.at[i, 'mcintosh_d'] = sk.mcintosh_d(counts)
    langdf.at[i, 'strong'] = sk.strong(counts)
    langdf.at[i, 'shannon'] = sk.shannon(counts, base=np.e)
    langdf.at[i, 'brillouin'] = sk.brillouin_d(counts)
    langdf.at[i, 'pielou'] = sk.pielou_e(counts)
    langdf.at[i, 'heip'] = sk.heip_e(counts)
    langdf.at[i, 'simpson_e'] = sk.simpson_e(counts)
    langdf.at[i, 'mcintosh_e'] = sk.mcintosh_e(counts)
    langdf.at[i, 'menhinick'] = sk.menhinick(counts)
    langdf.at[i, 'margalef'] = sk.margalef(counts)
    langdf.at[i, 'gini'] = sk.gini_index(counts)
    langdf.at[i, 'enspie'] = sk.enspie(counts)
    langdf.at[i, 'simpson'] = sk.simpson(counts)

# scale shannon
langdf['shannon_scaled'] = scale_minmax(langdf['shannon'])

# get time of day dataframe to grid
langrid = pd.merge(grid, langdf, how='outer', left_on='NRO', right_on='grid_id')
# ddrop empty rows
langrid = langrid.dropna(subset=['sents'])
# drop list column
langrid = langrid.drop(columns=['langs'])
# save to file
print('[INFO] - Saving results...')
langrid.to_file(args['output'], driver='GPKG')

# print done
print('[INFO] - ... done!')
