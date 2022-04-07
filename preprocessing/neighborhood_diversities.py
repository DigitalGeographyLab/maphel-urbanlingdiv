#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 12:10:25 2021

Calculate the linguistic diversities in grid cells for both Twitter and Instagram data.

@author: waeiski
"""
import geopandas as gpd
import pandas as pd
import pysal
import skbio.diversity.alpha as sk
import numpy as np

# read first language grid data in
grid = gpd.read_file('2015_diversities.gpkg')

# simplify data structure
grid = grid[['NRO', 'KUNTA', 'KUNTANRO', 'geometry']]

# convert grid id to integer to enable merging
grid['NRO'] = grid['NRO'].astype(int)

# define the neighbourhood criteria
gW = pysal.lib.weights.Queen.from_dataframe(grid, idVariable='NRO')

# filelists for Twitter and Instagram data
ifiles = ['insta2015_morning.pkl','insta2015_noon.pkl','insta2015_afternoon.pkl',
          'insta2015_evening.pkl','insta2015_night.pkl']
tfiles = ['twitter2015_morning.pkl','twitter2015_noon.pkl','twitter2015_afternoon.pkl',
          'twitter2015_evening.pkl','twitter2015_night.pkl']

# location lists of general points
ilocs = ['Helsinki, Finland', 'Espoo, Finland', 'Vantaa, Finland' ,'Helsinki',
         'Finland', 'Helsinki - Finland']

# read insta data in
for ix in range(len(ifiles)):
    print('[INFO] - Processing ' + str(ifiles[ix]) + '....')
    # get file and col
    file = ifiles[ix]
    # read timeofday df in
    df = pd.read_pickle('/filepath/to/socialmedia/' + file)
    # drop general pois
    df = df[~df['location_name'].isin(ilocs)]
    # join with grid
    df = gpd.sjoin(df, grid, op='intersects')
    # get grid ids
    idlist = df['NRO'].value_counts().index.tolist()
    # set language dictionary
    langdict = {}
    # loop over grid ids
    for gid in idlist:
        # get neighboring grid ids of current grid
        neigh_ids = list(gW[gid].keys())
        # extend list with current grid id
        neigh_ids.extend([gid])
        # fetch neighborhood grid cells for current grid id into dataframe
        ndf = df[df['NRO'].isin(neigh_ids)]
        # extract languages used in neighborhood grid cells
        langs = list(zip(ndf['language'].value_counts().index,
                         ndf['language'].value_counts().values))
        # add language use to language dictionary
        langdict[gid] = langs
    # generate a language use dataframe
    langdf = pd.DataFrame(langdict.items(), columns=['grid_id', 'langs'])
    # loop over language dataframe
    for i, row in langdf.iterrows():
        # get language counts from tuples
        counts = [x[1] for x in row['langs']]
        # calculate diversities
        langdf.at[i, 'sents'] = sum(counts)
        langdf.at[i, 'unique'] = sk.observed_otus(counts)
        langdf.at[i, 'singletons'] = sk.singles(counts)
        langdf.at[i, 'berger'] = sk.berger_parker_d(counts)
        langdf.at[i, 'dominance'] = sk.dominance(counts)
        langdf.at[i, 'mcintosh_d'] = sk.mcintosh_d(counts)
        langdf.at[i, 'strong'] = sk.strong(counts)
        langdf.at[i, 'shannon'] = sk.shannon(counts, base=np.e)
        langdf.at[i, 'brillouin'] = sk.brillouin_d(counts)
        langdf.at[i, 'pielou'] = sk.pielou_e(counts)
        langdf.at[i, 'heip'] = sk.heip_e(counts)
        langdf.at[i, 'simpson_e'] = sk.simpson_e(counts)
        langdf.at[i, 'mcintosh_e'] = sk.mcintosh_e(counts)
        langdf.at[i, 'menhinick'] = sk.menhinick(counts)
        langdf.at[i, 'margalef'] = sk.margalef(counts)
        langdf.at[i, 'gini'] = sk.gini_index(counts)
        langdf.at[i, 'enspie'] = sk.enspie(counts)
        
    # save current dataframe to pickle
    langdf.to_pickle('/results/' + file[:-4] + '_diversities.pkl')
    # get time of day dataframe to grid
    langrid = pd.merge(grid, langdf, how='outer', left_on='NRO', right_on='grid_id')
    # ddrop empty rows
    langrid = langrid.dropna(subset=['sents'])
    # drop list column
    langrid = langrid.drop(columns=['langs'])
    # save to file
    langrid.to_file('/results/diversities_' + file[:-4] + '.gpkg',driver='GPKG')

# read twitter in
for ix in range(len(tfiles)):
    print('[INFO] - Processing ' + str(tfiles[ix]) + '....')
    # get file and col
    file = tfiles[ix]
    # read timeofday df in
    df = pd.read_pickle('/filepath/to/socialmedia/' + file)
    # join with grid
    df = gpd.sjoin(df, grid, op='intersects')
    # drop bounding box coords
    df = df[df['locinfo_type'] == 'gps']
    # get grid ids
    idlist = df['NRO'].value_counts().index.tolist()
    # set language dictionary
    langdict = {}
    # loop over grid ids
    for gid in idlist:
        # get neighboring grid ids of current grid
        neigh_ids = list(gW[gid].keys())
        # extend list with current grid id
        neigh_ids.extend([gid])
        # fetch neighborhood grid cells for current grid id into dataframe
        ndf = df[df['NRO'].isin(neigh_ids)]
        # extract languages used in neighborhood grid cells
        langs = list(zip(ndf['language'].value_counts().index,
                         ndf['language'].value_counts().values))
        # add language use to language dictionary
        langdict[gid] = langs
    # generate a language use dataframe
    langdf = pd.DataFrame(langdict.items(), columns=['grid_id', 'langs'])
    # loop over language dataframe
    for i, row in langdf.iterrows():
        # get language counts from tuples
        counts = [x[1] for x in row['langs']]
        # calculate diversities
        langdf.at[i, 'sents'] = sum(counts)
        langdf.at[i, 'unique'] = int(sk.observed_otus(counts))
        langdf.at[i, 'singletons'] = sk.singles(counts)
        langdf.at[i, 'berger'] = sk.berger_parker_d(counts)
        langdf.at[i, 'dominance'] = sk.dominance(counts)
        langdf.at[i, 'mcintosh_d'] = sk.mcintosh_d(counts)
        langdf.at[i, 'strong'] = sk.strong(counts)
        langdf.at[i, 'shannon'] = sk.shannon(counts, base=np.e)
        langdf.at[i, 'brillouin'] = sk.brillouin_d(counts)
        langdf.at[i, 'pielou'] = sk.pielou_e(counts)
        langdf.at[i, 'heip'] = sk.heip_e(counts)
        langdf.at[i, 'simpson_e'] = sk.simpson_e(counts)
        langdf.at[i, 'mcintosh_e'] = sk.mcintosh_e(counts)
        langdf.at[i, 'menhinick'] = sk.menhinick(counts)
        langdf.at[i, 'margalef'] = sk.margalef(counts)
        langdf.at[i, 'gini'] = sk.gini_index(counts)
        langdf.at[i, 'enspie'] = sk.enspie(counts)
    # save current dataframe to pickle
    langdf.to_pickle('/results/' + file[:-4] + '_diversities.pkl')
    # get time of day dataframe to grid
    langrid = pd.merge(grid, langdf, how='outer', left_on='NRO', right_on='grid_id')
    # drop nones
    langrid = langrid.dropna(subset=['sents'])
    # drop list column
    langrid = langrid.drop(columns=['langs'])
    # save to file
    langrid.to_file('/results/diversities_' + file[:-4] + '.gpkg',driver='GPKG')
    
# print done
print('[INFO] - ... done!')