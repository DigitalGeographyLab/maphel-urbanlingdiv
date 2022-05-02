#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 15:54:08 2022

This script loads the pickled OLS models to create table of coefficients and p-values.

@author: waeiski
"""

import statsmodels.api as sm
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

# path to models
mdls = glob.glob(args['inputfolder'])

# empty list for coefficient dataframes
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

    # add to dataframe list
    mtab.append(comb)

# concatenate to construct one dataframe
result = pd.concat(mtab, axis=1)

# reorder columns
result = result[['morning', 'mo_pval', 'noon', 'no_pval', 'afternoon', 'af_pval', 'evening', 'ev_pval', 'night', 'ni_pval', 'fulldata', 'fu_pval']]

# print result in latex format for easier integration into latex
print(result.style.to_latex())
