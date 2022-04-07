#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 11:32:51 2021

@author: tuomvais
"""

import geopandas as gpd
from collections import Counter
import pandas as pd

# function to count dominant language of user
def count_dom(langlist):
    
    # get count of languages used
    total_use = len(langlist)
    
    # count unique languages
    uniques = len(set(langlist))
    
    # convert list to counted dictionary
    counts = Counter(langlist)
    
    # get most common language count
    lang_count = counts.most_common()[0][1]
    lang_code = counts.most_common()[0][0]
    
    # is multilingual
    if uniques == 1:
        langs = 'Monolingual'
    elif uniques == 2:
        langs = 'Bilingual'
    elif uniques > 2:
        langs = 'Multilingual'
    
    # count proportion
    prop = round((lang_count / total_use), 4)
    
    return lang_code, prop, uniques, langs

# read social media data in
origdf = gpd.read_file('/home/tuomvais/GIS/maphel_thirdplace/maphel_thirdplace/some_combined/combined/twinsta_2015_langid.gpkg')

# all users
allusers = origdf.groupby(by=['userid'])['language'].apply(list)

# group by user_id and list languages per user
df = origdf.groupby(by=['userid'])['language'].apply(list).reset_index(drop=False)

# drop users with only one language detection
df = df[df['language'].apply(lambda x: len(x) >= 5)]

# calculate dominant language properties
for i, row in df.iterrows():
    df.at[i, 'dom_lang'] = count_dom(row['language'])[0]
    df.at[i, 'lang_prop'] = count_dom(row['language'])[1]
    df.at[i, 'ulang_count'] = count_dom(row['language'])[2]
    df.at[i, 'monomulti'] = count_dom(row['language'])[3]

# get more clear cases of dominant language use
clear = df[df['lang_prop'] > 0.5]
domcounts = clear['dom_lang'].value_counts()

# get sentence counts based on filtered userids
sentdf = origdf.merge(clear, on='userid', how='inner')
sentcount = sentdf['language_x'].value_counts()

# get monolingual users
mono = df[df['monomulti'] == 'Monolingual']

# reproduce top language counts in the registry data
# registry values based on individuals
# some values based on sentences
from scipy.stats import spearmanr
spear = pd.DataFrame({'lang':'fi sv ru et so en ar zh ku sq fr ja ko es tr de'.split(' '),
                    'reg': [879011, 63903, 28404, 23169, 11735, 8616, 7562, 5753, 4830, 4252, 2052, 690, 384, 3353, 2777, 2318],
                    'some': [417531, 11242, 32375, 1174, 0, 233881, 944, 410, 0, 0, 9761, 5961, 2874, 2736, 2453, 2160]})

# registry values based on individuals
# some values based on users from domcounts data from above
spear = pd.DataFrame({'lang':'fi sv ru et so en ar zh ku sq fr ja ko es tr de'.split(' '),
                    'reg': [879011, 63903, 28404, 23169, 11735, 8616, 7562, 5753, 4830, 4252, 2052, 690, 384, 3353, 2777, 2318],
                    'some': [29742, 1171, 5025, 78, 0, 18614, 88, 91, 0, 0, 93, 750, 459, 279, 96, 122]})

# calculate spearman rank
rho, p = spearmanr(spear['reg'], spear['some'])

# print results
print('Distribution of mono/multilingual users: \n')
print(clear['monomulti'].value_counts())
print('\n')
print('Top 10 main languages of users: \n')
print(clear['dom_lang'].value_counts()[:10])
print('\n')
print('Top 10 languages of monolingual users: \n')
print(mono['dom_lang'].value_counts()[:10])
