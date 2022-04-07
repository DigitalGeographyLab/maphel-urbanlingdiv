#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 15:54:08 2022

This script loads the pickled OLS models to print model summaries.

@author: waeiski
"""
import statsmodels.api as sm
import pandas as pd
import geopandas as gpd
import glob

# path to models
path = '/home/waeiski/GIS/maphel_thirdplace/ols/*.pkl'

# list all pickled models
mdls = glob.glob(path)

# read models
mo = sm.load(mdls[1])
no = sm.load(mdls[3])
af = sm.load(mdls[0])
ev = sm.load(mdls[4])
ni = sm.load(mdls[2])
fu = sm.load(mdls[5])

# print model summaries as latex
print(mo.summary().as_latex())

print(no.summary().as_latex())

print(af.summary().as_latex())

print(ev.summary().as_latex())

print(ni.summary().as_latex())

print(fu.summary().as_latex())
