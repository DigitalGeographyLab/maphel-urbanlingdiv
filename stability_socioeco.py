#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 09:57:23 2021

@author: tuomvais
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np

# read in rtk data
rtk = gpd.read_file('/home/waeiski/GIS/maphel_thirdplace/RTK_data_4_posts_users.gpkg')

# read in stability data of low and high diversity clusters
histab = gpd.read_file('/home/waeiski/GIS/maphel_thirdplace/stability_high.gpgk', driver='GPKG')
lowstab = gpd.read_file('/home/waeiski/GIS/maphel_thirdplace/stability_low.gpgk', driver='GPKG')

# classify stabilities for low stability data
lowstab['stabclass'] = lowstab['stability'].apply(lambda x: 'High' if x >= 4 else ('Moderate' if x == 3 else 'Low'))

# add classifications for stability data
histab['divclass'] = 'High diversity'
lowstab['divclass'] = 'Low diversity'

# join stability data together
stab = histab.append(lowstab).reset_index(drop=True)

# spatial join stability data with rtk data
joined = stab.sjoin(rtk, how='inner', predicate='within')

# set attributes
attrs = ['he_vakiy', 'he_naiset', 'he_miehet', 'he_kika', 'kids_u18',
         'pensioners', 'mean_income', 'low_income_prop', 'med_income_prop',
         'hi_income_prop', 'employed_prop', 'unemployed_prop', 'students_prop',
         'basic_ed_prop', 'high_school_prop', 'vocational_ed_prop',
         'bachelor_ed_prop', 'masters_ed_prop', 'med_ed_prop', 'high_ed_prop',
         'total_jobs', 'mean_livingspace', 'avg_sq_m', 'res_buildings_prop',
         'oth_buildings_prop']


# set y labels
labels = ['Total population', 'Total women', 'Total men', 'Average age', 'Kids under 18 years',
        'Pensioners', 'Mean income (€)', 'Low income individuals (%)', 'Medium income individuals (%)', 'High income individuals (%)',
        'Employed (%)', 'Unemployed (%)', 'Students (%)', 'Basic education level (%)', 'High school education level (%)',
        'Vocational education level (%)', 'BSc education level (%)', 'MSc education level (%)',
        'Medium education level (%)', 'High education level (%)', 'Total jobs', 'Mean livingspace (m\u00b2)',
        'Average m\u00b2', 'Residential buildings (%)', 'Non-residential buildings (%)']

# initiate plotting
sns.set(font_scale=1.5)
colors = sns.color_palette(['#4dac26', '#d8a600', '#1b21d0'])

# loop over attributes to plot
for i, attr in enumerate(attrs):
    
    # create figure
    plt.figure(figsize=(8,14))
    
    # plot categorical plot
    g = sns.catplot(x='stabclass', y=attr, hue='stabclass', kind='boxen',
                    dodge=False, col='divclass', palette=colors,
                    order=['High', 'Moderate', 'Low'],
                    hue_order=['High', 'Moderate', 'Low'], data=joined)
    
    # final details
    (g.set_axis_labels("", "")
     .set_xticklabels(['High', 'Moderate', 'Low'])
     .set_titles("{col_name}").despine(left=True))
    
    # get flexible text location value
    tloc = - joined[attr].max() / 5.2
    
    # add text
    plt.text(-4.2, tloc, "Temporal\nstability:", fontsize=12)
        
    # save plot
    plt.savefig('/home/tuomvais/GIS/maphel_thirdplace/plots/div_stab_' + attr + '_boxen.pdf' , dpi=300, bbox_inches='tight')


# INCREASE FONT SIZE PLOTS

# set attributes
attrs = ['low_income_prop', 'med_income_prop',
         'hi_income_prop', 'employed_prop', 'unemployed_prop', 'students_prop',
         'basic_ed_prop', 'high_school_prop', 'vocational_ed_prop',
         'bachelor_ed_prop', 'masters_ed_prop', 'med_ed_prop', 'high_ed_prop',
         'res_buildings_prop', 'oth_buildings_prop', 'small_house_prop', 'block_house_prop',
         'owner_household_prop', 'rent_household_prop', 'horeca_jobs', 'construction_jobs']

# initiate plotting
sns.set(font_scale=1.7)
colors = sns.color_palette(['#4dac26', '#d8a600', '#1b21d0'])

# loop over attributes to plot
for i, attr in enumerate(attrs):
    
    # create figure
    plt.figure(figsize=(8,14))
    
    # plot categorical plot
    g = sns.catplot(x='stabclass', y=attr, hue='stabclass', kind='boxen',
                    dodge=False, col='divclass', palette=colors,
                    order=['High', 'Moderate', 'Low'],
                    hue_order=['High', 'Moderate', 'Low'], data=joined)
    
    # final details
    (g.set_axis_labels("", "")
     .set_xticklabels(['High', 'Moderate', 'Low'])
     .set_titles("{col_name}").despine(left=True).set(ylim=(-0.05,1.05)))
    
    # add text
    plt.text(-4.2, -0.18,u"Temporal\nstability:", fontsize=12)
    
    # save plot
    plt.savefig('/home/waeiski/GIS/maphel_thirdplace/plots/div_stab_' + attr + '_boxen.pdf' , dpi=300, bbox_inches='tight')


# attributes identified with OLS/GWR
attrs = ['he_kika', 'unemployed_prop', 'masters_ed_prop', 'block_house_prop',
         'rent_household_prop', 'vocational_ed_prop', 'low_income_prop',
         'horeca_jobs', 'construction_jobs']

# initiate plotting
sns.set(font_scale=1.7)
colors = sns.color_palette(['#4dac26', '#d8a600', '#1b21d0'])

# loop over attributes to plot
for i, attr in enumerate(attrs):
    
    # create figure
    plt.figure(figsize=(8,14))
    
    # plot categorical plot
    g = sns.catplot(x='stabclass', y=attr, hue='stabclass', kind='box',
                    dodge=False, col='divclass', palette=colors,
                    order=['High', 'Moderate', 'Low'],
                    hue_order=['High', 'Moderate', 'Low'], data=joined)
    
    g = sns.catplot(x='stabclass', y=attr, kind='swarm', palette=colors,
                    col='divclass', order=['High', 'Moderate', 'Low'],
                    hue_order=['High', 'Moderate', 'Low'],
                    dodge=False, data=joined)
    
    # final details
    (g.set_axis_labels("", "")
     .set_xticklabels(['High', 'Moderate', 'Low'])
     .set_titles("{col_name}").despine(left=True).set(ylim=(-0.05,1.05)))
    
    # add text
    plt.text(-4.2, -0.18,u"Temporal\nstability:", fontsize=12)
    
    # save plot
    plt.savefig('/home/waeiski/GIS/maphel_thirdplace/plots/div_stab_' + attr + '_boxen.pdf' , dpi=300, bbox_inches='tight')




















# set attributes
attrs = ['he_vakiy', 'he_naiset', 'he_miehet', 'he_kika', 'kids_u18',
         'pensioners', 'mean_income',
         'total_jobs', 'mean_livingspace', 'avg_sq_m']

# initiate plotting
sns.set(font_scale=1.7)
colors = sns.color_palette(['#4dac26', '#d8a600', '#1b21d0'])

# loop over attributes to plot
for i, attr in enumerate(attrs):
    
    # create figure
    plt.figure(figsize=(8,14))
    
    # plot categorical plot
    g = sns.catplot(x='stabclass', y=attr, hue='stabclass', kind='boxen',
                    dodge=False, col='divclass', palette=colors,
                    order=['High', 'Moderate', 'Low'],
                    hue_order=['High', 'Moderate', 'Low'], data=joined)
    
    # final details
    (g.set_axis_labels("", "")
     .set_xticklabels(['High', 'Moderate', 'Low'])
     .set_titles("{col_name}").despine(left=True))
    
    # get flexible text location value
    tloc = - joined[attr].max() / 5.2
    
    # add text
    plt.text(-4.2, tloc, "Temporal\nstability:", fontsize=12)
    
    # save plot
    plt.savefig('/home/waeiski/GIS/maphel_thirdplace/plots/div_stab_' + attr + '_boxen.pdf' , dpi=300, bbox_inches='tight')



##################################################################################
# old stuff
# loop over enumerated attributes
for i, atr in enumerate(attrs):
    
    # plot distributions
    sns.catplot(x='stabclass', y='vocational_ed_prop', hue='stabclass', kind='boxen',
                dodge=False, col='divclass', order=['High', 'Moderate', 'Low'],
                data=joined, scale='linear')
    
    # set titles according to alphabets
    axes[i].set_title(titles[i], loc='left', fontsize=16, y=0.9, x=0.03)
    
    # clarify labeling of y axis
    axes[i].set_ylabel(labs[i])
    
    # drop x labels
    axes[i].set_xlabel('')
    
    



# create empty list for file paths
files = []

# populate list with geopackage file paths
for gpkg in glob.glob('/home/tuomvais/GIS/maphel_thirdplace/socioeco_analysis/*.gpkg'):
    if 'RTK' in gpkg:
        files.append(gpkg)
    else:
        pass

# read file in
df = gpd.read_file('/home/tuomvais/GIS/maphel_thirdplace/socioeco_analysis/joined_RTK_stability_clusters2.gpkg')

# drop rows without clusters
df = df.dropna(subset=['stabclass'])

# get only the separate cluster files
files = files[:-1]

# list of attributes to plot
attrs = ['hi_income_prop', 'med_income_prop', 'low_income_prop', 'unemployed_prop',
         'high_ed_prop', 'med_ed_prop', 'basic_ed_prop', 'mean_livingspace']
attrs = ['hi_income_prop', 'med_income_prop', 'low_income_prop', 'students_prop',
         'res_buildings_prop', 'oth_buildings_prop', 'total_jobs', 'test']

# list of y labels
labs = ['High income (%)', 'Middle income (%)', 'Low income (%)', 'Unemployed (%)',
        'University education (%)', 'High school education (%)', 'Basic education (%)',
        'Mean livingspace (m\u00b2)']
labs = ['High income (%)', 'Middle income (%)', 'Low income (%)', 'Students (%)',
        'University education (%)', 'High school education (%)', 'Basic education (%)',
        'Mean income']

# list of titles
titles = ['a.','b.','c.','d.','e.','f.','g.','h.']

# test standardizing square meters
sqm = (df['avg_sq_m'] - df['avg_sq_m'].mean()) / df['avg_sq_m'].std()
df['test'] = sqm

# initiate plotting
sns.set()
fig, axes = plt.subplots(2, 4, figsize=(20,10))
axes = axes.flatten()

# loop over enumerated attributes
for i, atr in enumerate(attrs):
    
    # plot distributions
    sns.boxplot(ax=axes[i], x='stabclass', y=atr,  dodge=False,
                   order=['High', 'Moderate', 'Low'], data=df)
    
    # set titles according to alphabets
    axes[i].set_title(titles[i], loc='left', fontsize=16, y=0.9, x=0.03)
    
    # clarify labeling of y axis
    axes[i].set_ylabel(labs[i])
    
    # drop x labels
    axes[i].set_xlabel('')

plt.savefig('/home/tuomvais/GIS/maphel_thirdplace/plots/stability_socioeco.pdf', dpi=300, bbox_inches='tight')
    