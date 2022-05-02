#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 12:16:36 2021

This script runs the local bivariate Moran's I analysis with Shannon entropy and Simpson diversity.

@author: waeiski
"""
from pysal.lib import weights
from pysal.explore import esda
import geopandas as gpd
import pandas as pd
import glob
from sklearn.preprocessing import MinMaxScaler
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Get path to input folder
ap.add_argument("-if", "--inputfolder", required=True,
                help="Path to folder containing grid files with diversity values.")

# Get path to input folder
ap.add_argument("-of", "--outputfolder", required=True,
                help="Path to output folder.")

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

# diversity indices
divs = ['shannon', 'shannon_scaled', 'margalef', 'unique', 'simpson']

# create empty list for file paths
files = []

# populate list with geopackage file paths
for gpkg in glob.glob(args['inputfolder'] + '*.gpkg'):
    files.append(gpkg)

# loop over the files
print('[INFO] - Starting processing...')
for file in files:

    # get file name
    fn = file.split('/')[-1]

    # read file in
    print('[INFO] - Processing ' + fn + '...')
    grid = gpd.read_file(file)

    # convert grid id to integer
    grid['NRO'] = grid['NRO'].astype(int)

    # get contiguity weights
    #gW = weights.Queen.from_dataframe(grid, idVariable='NRO', silence_warnings=True)
    #gW = weights.Kernel.from_dataframe(grid, ids='NRO', fixed=False, k=8, silence_warnings=True)
    #gW = weights.DistanceBand.from_dataframe(grid, ids='NRO', threshold=400, binary=False, silence_warnings=True)
    gW = weights.KNN.from_dataframe(grid, ids='NRO', k=8, silence_warnings=True)

    # try to unify naming scheme
    try:
        grid['unique'] = grid['unique_langs']
    except:
        pass

    # calculate simpson diversity
    grid['simpson'] = 1 - grid['dominance']

    # scale shannon entropy to 0,1 scale
    grid['shannon_scaled'] = scale_minmax(grid['shannon'])

    # loop over diversities
    for div in divs:

        # get output column name
        outcol = div + '_cl'

        #  calculate local moran's i
        li = esda.Moran_Local(grid[div], gW, transformation='r', permutations=9999, n_jobs=3)

        # how many are statistically significant
        sigs = (li.p_sim < 0.01).sum()

        # get all signifiant cells
        sig = 1 * (li.p_sim < 0.001)

        # get hot spots (HH)
        hotspot = 1 * (sig * li.q==1)

        # get cold spots (LL)
        coldspot = 3 * (sig * li.q==3)

        # get cold outliers (LH)
        lohi = 2 * (sig * li.q==2)

        # get hot outliers (HL)
        hilo = 4 * (sig * li.q==4)

        # combine statistically significant spots
        spots = hotspot + coldspot + lohi + hilo

        # add column and values
        grid[outcol] = spots

    # get output column name
    outBVcol = 'sha_sim_cl'

    #  calculate local moran's i
    li = esda.Moran_Local_BV(grid['shannon_scaled'], grid['simpson'], gW, transformation='r', permutations=9999)

    # how many are statistically significant
    sigs = (li.p_sim < 0.01).sum()

    # get all signifiant cells
    sig = 1 * (li.p_sim < 0.001)

    # get hot spots (HH)
    hotspot = 1 * (sig * li.q==1)

    # get cold spots (LL)
    coldspot = 3 * (sig * li.q==3)

    # get cold outliers (LH)
    lohi = 2 * (sig * li.q==2)

    # get hot outliers (HL)
    hilo = 4 * (sig * li.q==4)

    # combine statistically significant spots
    spots = hotspot + coldspot + lohi + hilo

    # add column and values
    grid[outBVcol] = spots

    # save to file
    grid.to_file(args['output'] + 'knn8_cl_' + fn, driver='GPKG')

# print message
print('[INFO] - ... done!')
