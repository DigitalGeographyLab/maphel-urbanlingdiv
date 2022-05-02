#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 16:12:48 2022

This script reads the VIF test results and prints out the scores.

@author: waeiski
"""

import pandas as pd
import geopandas as gpd
import glob
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Get path to input folder
ap.add_argument("-if", "--inputfolder", required=True,
                help="Path to folder containing pickled model files.")

# parse arguments
args = vars(ap.parse_args())

# path to times of day data
path = args['inputfolder']

# list vif files
files = glob.glob(path + '*.pkl')

# empty lists for initial and final files
invif = []
finvif = []

# loop over files
for file in files:
    # check if initial or final
    if 'initial' in file:
        invif.append(file)
    elif 'final' in file:
        finvif.append(file)

# loop over initial files
for file in invif:
    # get time of day
    timeofday = file.split('_')[-2]
    # read file
    df = pd.read_pickle(file)
    # print contents
    print('[INFO] - Now printing initial VIF scores for ' + str(timeofday))
    print(df)

# loop over final files
for file in finvif:
    # get time of day
    timeofday = file.split('_')[-2]
    # read file
    df = pd.read_pickle(file)
    # print contents
    print('[INFO] - Now printing final VIF scores for ' + str(timeofday))
    print(df)

print('[INFO] - ... done!')
