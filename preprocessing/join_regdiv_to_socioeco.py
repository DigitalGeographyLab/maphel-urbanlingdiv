#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 14:38:08 2022

@author: tuomvais
"""

import geopandas as gpd
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Get path to register input file
ap.add_argument("-r", "--register", required=True,
                help="Path to input register file (type geopackage).")

# Get path to grid database input file
ap.add_argument("-i", "--input", required=True,
                help="Path to input RTK database file (type geopackage).")

# Get path to output file
ap.add_argument("-o", "--output", required=True,
                help="Path to output file (type geopackage).")

# parse arguments
args = vars(ap.parse_args())

# read register data with diversities in
reg = gpd.read_file(args['register'])

# read the socio-economic data in
rtk = gpd.read_file(args['input'])

# simplify data structure in register data
reg = reg[['NRO','shannon']]

# rename shannon
reg = reg.rename(columns={'shannon':'reg_shannon'})

# join data
joined = rtk.merge(reg, on='NRO')

# save to disk
joined.to_file(args['output'], driver='GPKG')
