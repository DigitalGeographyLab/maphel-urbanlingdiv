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
rtk = gpd.read_file('RTK_data.gpkg')

# read in stability data of low and high diversity clusters
histab = gpd.read_file('stability_high.gpgk', driver='GPKG')
lowstab = gpd.read_file('stability_low.gpgk', driver='GPKG')

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
    plt.savefig('plots/div_stab_' + attr + '_boxen.pdf' , dpi=300, bbox_inches='tight')