#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 14:42:07 2022

OLS regression analysis for all times of day with all cells, not just with
identified cluster locations.

@author: Tuomas Väisänen
"""

import numpy as np
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy.stats import spearmanr
import geopandas as gpd
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import random

# function to indicate significant correlations between variables
def corr_vars(dataframe, x_col, var_cols, threshold=0.05):
    
    # fill na
    dataframe = dataframe.fillna(value=0)
    
    # create df list for significant vars
    variables = pd.DataFrame(columns=['variable', 'corr', 'pvalue'])
    
    # loop over columns
    for i in range(len((var_cols))):
        cor, p = spearmanr(dataframe[x_col].values, dataframe[var_cols[i]].values)
        
        # add to df
        variables.at[i, 'variable'] = var_cols[i]
        variables.at[i, 'corr'] = cor
        variables.at[i, 'pvalue'] = p
        
    # order dataframe by p-value
    variables = variables.sort_values(by='corr').reset_index(drop=True)
    
    # return result
    return variables
    

# backward feature selection function for regression
def backward_regression(X, y,
                           initial_list=[], 
                           threshold_in=0.01, 
                           threshold_out = 0.05, 
                           verbose=True):
    # include all dataframe columns for feature selection
    included = list(X.columns)
    
    # while loop until p-value thresholds are met
    while True:
        changed = False
        
        # create an ols regression model and fit it
        model = sm.OLS(y, sm.add_constant(pd.DataFrame(X[included])), missing='drop').fit()
        
        # use all coefs except intercept
        pvalues = model.pvalues.iloc[1:]
        
        # extract worst p value
        worst_pval = pvalues.max() # null if pvalues is empty
        
        # check if worst pvalue exceeds threshold
        if worst_pval > threshold_out:
            
            # indicate change in feature list
            changed = True
            
            # select worst feature
            worst_feature = pvalues.idxmax()
            
            # drop worst feature from list
            included.remove(worst_feature)
            
            # indicate which feature was dropped and why
            if verbose:
                print('Drop {} with p-value {}'.format(worst_feature, worst_pval))
                
        # if no features were dropped, feature selection is complete and return list
        if not changed:
            break
        
    return included

# function for variance inflation factor filtering
def vif_filter(dataframe, cols, initial_list=[], threshold=5):
    
    # include all columns
    included = list(dataframe[cols].columns)
    
    # while loop until vif thresholds are met
    while True:
        changed = False
        
        # get independent variables and drop nans
        X = dataframe[included].dropna()
        
        # create vif dataframe
        vif = pd.DataFrame()
        
        # add feature column
        vif['feature'] = X.columns
        
        # calculate vif
        vif['VIF'] = [variance_inflation_factor(X.values, i) for i in range(len(X.columns))]
        
        # get worst vif
        worst_vif = vif['VIF'].max()
        
        # check against threshold
        if worst_vif > threshold:
            
            # indicate change
            changed = True
            
            # get name of worst feature
            worst = vif.loc[vif['VIF'].argmax()]['feature']
            
            # drop worst feature from included
            included.remove(worst)
            print('[INFO] - Drop {} with VIF score {}'.format(worst, worst_vif))
        
        # break loop if nothing changed
        if not changed:
            
            # for final vif dataframe get independent variables and drop nans
            X = dataframe[included].dropna()
            
            # create vif dataframe
            vif = pd.DataFrame()
            
            # add feature column
            vif['feature'] = X.columns
            
            # calculate vif
            vif['VIF'] = [variance_inflation_factor(X.values, i) for i in range(len(X.columns))]
            
            break
    
    # return trimmed list
    return included, vif

# read in grid database (rtk) data
rtk = gpd.read_file('/RTK_data.gpkg', driver='GPKG')
rtk = rtk.set_crs('epsg:3067', allow_override=True)

# combine job data
rtk['horeca_retail'] = rtk['horeca_jobs'] + rtk['retail_jobs']
rtk['postcount_scaled'] = rtk['postcount'] / rtk['postcount'].max()
rtk['bluecollar'] = rtk['basic_ed_prop'] + rtk['vocational_ed_prop']
rtk['reg_shan_scaled'] = rtk['reg_shannon'] / rtk['reg_shannon'].max()
rtk['total_jobs_scaled'] = rtk['total_jobs'] / rtk['total_jobs'].max()

# read in full twinsta data
twinsta = gpd.read_file('/twinsta_diversities_2015_no_timeofday.gpkg', driver='GPKG')

# join full data with rtk data
df = twinsta.sjoin(rtk, how='inner', predicate='within')

# dependent variables
y = df['shannon_scaled']

# these vars decided with olle and used on 28.01 (r2 0.179)
vcols =['he_kika', 'hi_income_prop', 'unemployed_prop', 'students_prop',
        'basic_ed_prop','vocational_ed_prop', 'high_ed_prop', 'oth_buildings_prop',
        'block_dwel_prop', 'rent_household_prop', 'edsci_jobs',
        'entertainment_jobs', 'horeca_retail', 'intorg_jobs', 'total_jobs_scaled',
        'postcount_scaled', 'reg_shan_scaled']

# ensure variable order is not affecting results
random.shuffle(vcols)

# run vif filter
vif_vars, vif_df = vif_filter(df, vcols, threshold=5)
vif_df.to_pickle('/ols/vif/initial_fulldata_vif.pkl')

# randomize variable order to ensure order is not affecting the analysis
random.shuffle(vif_vars)

# run backwards regression
reg_vars = backward_regression(df[vif_vars], y, threshold_out=0.1)

# get final vif scores
vif_list, vifs = vif_filter(df, reg_vars, threshold=5)
vifs.to_pickle('/ols/vif/final_fulldata_vif.pkl')

# add constants to statsmodel
X_var = sm.add_constant(df[reg_vars[:]])

# fit regression
model = sm.OLS(y, X_var, missing='drop',).fit()

# plot residuals
plt.figure(figsize=(11,9))
m = sns.scatterplot(x=model.fittedvalues, y=model.resid, legend=True)
m.set_title('OLS residuals')
m.set_xlabel('Fitted value')
m.set_ylabel('Residual')
plt.savefig('ols/plots/fulldata_residuals.png')

# print summary
model.summary()

# save model
model.save('/ols/ols_fulldata.pkl')

# save dataframe
df.to_file('/ols/geopackages/full_data.gpkg', driver='GPKG')


# print latex formatted summary
print(model.summary().as_latex())


##############################################################################
##############################################################################
# Here's the time of day version

# path to times of day data
path = '/autocorr/'

# get list of times-of-day data
timesofday = glob.glob(path + '/knn8*.gpkg')

# filter list to contain only social media data
timesofday = [i for i in timesofday if not ('knn8_cl_2015' in i)]

# loop over filenames
for file in timesofday:
    
    # extract time
    timed = file.split('_')[-2]
    
    # read file in
    print('[INFO] - Processing data from {}...\n'.format(timed))
    timedf = gpd.read_file(file)
    
    # spatial join with socioeco data
    timedf = timedf.sjoin(rtk, how='inner', predicate='within')
    timedf.to_pickle('/ols/dataframes/{}_df.pkl'.format(timed))
    timedf.to_file('/ols/geopackages/{}_df.gpkg'.format(timed), driver='GPKG')
    
    # get dependent variable
    y = timedf['shannon_scaled']
    
    # copy original variables as separate list
    dynvars = vcols.copy()
    
    # get correct dynamic pop field name
    dynpop = 'scaled_' + timed
    
    # add dynpop field to variable list
    dynvars.append(dynpop)
    
    # randomize variable order
    random.shuffle(dynvars)
    
    # run vif filter
    vif_vars, vifdf = vif_filter(timedf, dynvars, threshold=5)
    vifdf.to_pickle('/ols/vif/initial_{}_vif.pkl'.format(timed))
    
    # randomize vif variable order
    random.shuffle(vif_vars)
    
    # run backwards regression on filtered variables
    reg_vars = backward_regression(timedf[vif_vars], y)
    
    # get final vif scores and save to df
    viflist, vifs = vif_filter(timedf, reg_vars, threshold=5)
    vifs.to_pickle('/ols/vif/final_{}_vif.pkl'.format(timed))
    
    # add constant to model
    X_var = sm.add_constant(timedf[reg_vars])
    
    # fit regression on time of day data
    dynmodel = sm.OLS(y, X_var, missing='drop').fit()
    
    # plot residuals
    plt.figure(figsize=(11,9))
    m = sns.scatterplot(x=dynmodel.fittedvalues, y=dynmodel.resid, legend=True)
    m.set_title('OLS residuals during {}'.format(timed))
    m.set_xlabel('Fitted value')
    m.set_ylabel('Residual')
    plt.savefig('/ols/plots/{}_residuals.png'.format(timed))
    
    # print summary
    print('\n[INFO] - For ' + timed + ' the R squared: ' + str(dynmodel.rsquared) + '\n')
    dynmodel.summary()
    
    # save model
    dynmodel.save('/ols/ols_{}.pkl'.format(timed))
    
    # print summary with latex formatting
    print(dynmodel.summary().as_latex())