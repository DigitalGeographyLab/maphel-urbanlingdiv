#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 14:38:08 2022

@author: tuomvais
"""

import geopandas as gpd

# read register data with diversities in
reg = gpd.read_file('2015_diversities.gpkg')

# read the socio-economic data in
rtk = gpd.read_file('RTK_data.gpkg')

# simplify data structure in register data
reg = reg[['NRO','shannon']]

# rename shannon
reg = reg.rename(columns={'shannon':'reg_shannon'})

# join data
joined = rtk.merge(reg, on='NRO')

# save to disk
joined.to_file('RTK_data.gpkg')