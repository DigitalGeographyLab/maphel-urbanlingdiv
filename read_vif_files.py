#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 16:12:48 2022

@author: waeiski
"""

import pandas as pd
import geopandas as gpd
import glob

# path to times of day data
path = '/home/waeiski/GIS/maphel_thirdplace/ols/vif/'

# get list of times-of-day data
invif = glob.glob(path + 'initial*.pkl')
finvif = glob.glob(path + 'final*.pkl')

# initial data
inmo = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/initial_morning_vif.pkl')
inno = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/initial_noon_vif.pkl')
inaf = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/initial_afternoon_vif.pkl')
inev = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/initial_evening_vif.pkl')
inni = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/initial_night_vif.pkl')

# final data
fimo = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/final_morning_vif.pkl')
fino = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/final_noon_vif.pkl')
fiaf = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/final_afternoon_vif.pkl')
fiev = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/final_evening_vif.pkl')
fini = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/final_night_vif.pkl')