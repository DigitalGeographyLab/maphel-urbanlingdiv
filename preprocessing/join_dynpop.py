#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script joins the dynamic population data to the grid database.

@author: tuomvais
"""

import pandas as pd
import geopandas as gpd
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Get path to input csv file
ap.add_argument("-c", "--csv", required=True,
                help="Path to input csv file.")

ap.add_argument("-i", "--input", required=True,
                help="Path to input dynamic population file (type geojson).")

# Get path to output file
ap.add_argument("-o", "--output", required=True,
                help="Path to output file (type geopackage).")

# parse arguments
args = vars(ap.parse_args())

# read dynamic population data in
df = pd.read_csv(args['csv'])

# read target grid in
grid = gpd.read_file(args['input'])

# summarize pop across times of day
df['pop_morning'] = df[['H6','H7','H8','H9']].sum(axis=1)
df['pop_noon'] = df[['H10','H11','H12','H13']].sum(axis=1)
df['pop_afternoon'] = df[['H14','H15','H16','H17']].sum(axis=1)
df['pop_evening'] = df[['H18','H19','H20','H21']].sum(axis=1)
df['pop_night'] = df[['H22','H23','H0','H1','H2','H3','H4','H5']].sum(axis=1)

# scale according to population maximum
df['morning_scaled'] = df['pop_morning'] / df['pop_morning'].max()
df['noon_scaled'] = df['pop_noon'] / df['pop_noon'].max()
df['afternoon_scaled'] = df['pop_afternoon'] / df['pop_afternoon'].max()
df['evening_scaled'] = df['pop_evening'] / df['pop_evening'].max()
df['night_scaled'] = df['pop_night'] / df['pop_night'].max()

# join to grid
joined = grid.merge(df, on='YKR_ID')

# simplify data structure
joined = joined[['pop_morning', 'pop_noon', 'pop_afternoon', 'pop_evening',
                 'pop_night','YKR_ID', 'geometry']]

# get centroids
joined['centroid'] = joined.centroid

# drop geometry and rename centroid to geometry
joined = joined.drop(columns=['geometry'])
joined = joined.rename(columns={'centroid':'geometry'})

# read grid database in
rtk = gpd.read_file('RTK_data.gpkg')

# conduct spatial join
final = joined.sjoin(rtk, how="right")

# save to file
final.to_file(args['output'], driver='GPKG')
