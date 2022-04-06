#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 15:54:08 2022

@author: waeiski
"""

import statsmodels.api as sm
import pandas as pd
import geopandas as gpd
import glob

# path to models
path = '/home/waeiski/GIS/maphel_thirdplace/ols/*.pkl'

mdls = glob.glob(path)

# empty dataframe for coefficients
mtab = []

# loop over models
for model in mdls:
    
    # extract indicator for model
    indi = model.split('_')[-1][:-4]
    
    # read model
    m = sm.load(model)
    
    print('The adjusted R squared for ' + indi + ' is: ' + str(round(m.rsquared_adj, 4)))
    
    # get coefficients
    coef = round(m.params, 3)
    
    # get pvalues
    pval = round(m.pvalues, 4)
    
    # convert to dataframes
    coef = coef.rename(indi).reset_index()
    pval = pval.rename(indi[:2] + '_pval').reset_index()
    
    # harmonize dynamic population naming scheme
    coef['index'] = coef['index'].str.replace('scaled_' + indi, 'dyn_pop')
    pval['index'] = pval['index'].str.replace('scaled_' + indi, 'dyn_pop')
    
    # combine
    comb = coef.merge(pval, on='index').set_index('index')
    
    # add to datafrmae
    mtab.append(comb)
    
# construct dataframes
result = pd.concat(mtab, axis=1)

# reorder columns
result = result[['morning', 'mo_pval', 'noon', 'no_pval', 'afternoon', 'af_pval', 'evening', 'ev_pval', 'night', 'ni_pval', 'fulldata', 'fu_pval']]

# print result in latex format
print(result.style.to_latex())

# read models
mo = sm.load(mdls[1])
no = sm.load(mdls[3])
af = sm.load(mdls[0])
ev = sm.load(mdls[4])
ni = sm.load(mdls[2])
fu = sm.load(mdls[5])

print(mo.summary().as_latex())

print(no.summary().as_latex())

print(af.summary().as_latex())

print(ev.summary().as_latex())

print(ni.summary().as_latex())

print(fu.summary().as_latex())
