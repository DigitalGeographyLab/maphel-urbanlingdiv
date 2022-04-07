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

# path to times of day data
path = '/home/waeiski/GIS/maphel_thirdplace/ols/vif/'

# get list of times-of-day data
invif = glob.glob(path + 'initial*.pkl')
finvif = glob.glob(path + 'final*.pkl')

# initial data
inmo = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/initial_morning_vif.pkl')
print(inmo)
inno = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/initial_noon_vif.pkl')
print(inno)
inaf = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/initial_afternoon_vif.pkl')
print(inaf)
inev = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/initial_evening_vif.pkl')
print(inev)
inni = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/initial_night_vif.pkl')
print(inni)

# final data
fimo = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/final_morning_vif.pkl')
print(fimo)
fino = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/final_noon_vif.pkl')
print(fino)
fiaf = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/final_afternoon_vif.pkl')
print(fiaf)
fiev = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/final_evening_vif.pkl')
print(fiev)
fini = pd.read_pickle('/home/waeiski/GIS/maphel_thirdplace/ols/vif/final_night_vif.pkl')
print(fini)