#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:16:42 2021

This script combines Twitter and Instagram data from 2015

@author: waeiski
"""
import geopandas as gpd

# read social media data in
twitter = gpd.read_file('/twitter2015_langid.gpkg')
instagram = gpd.read_file('insta2015_langid_tm35fin.gpkg')

# generic insta locations to drop to reduce post accumulation to single geographical point
ilocs = ['Helsinki, Finland', 'Espoo, Finland', 'Vantaa, Finland' ,'Helsinki',
         'Finland', 'Helsinki - Finland']

# drop the aforementioned locations
instagram = instagram[~instagram['location_name'].isin(ilocs)]

# simplify data structures
twitter = twitter[['id', 'userid', 'time_local', 'hour', 'weekday', 'sents', 'language', 'prob', 'charlen', 'geometry']]
instagram = instagram[['id', 'userid', 'time_local', 'hour', 'weekday', 'sents', 'language', 'prob', 'charlen', 'geometry']]

# rename 
twitter = twitter.rename(columns={'userid':'userid'})

# save platform info
twitter['platform'] = 'twitter'
instagram['platform'] = 'instagram'

# join twitter and instagram data
joined = twitter.append(instagram, ignore_index=True)

# get coordinates
joined['x_coord'] = joined['geometry'].x
joined['y_coord'] = joined['geometry'].y

# drop duplicate cross-posts at identical locations
joined = joined.drop_duplicates(subset=['sents','x_coord','y_coord'])

# separate data by times of day
morning = joined[joined['hour'].isin([6,7,8,9])]
noon = joined[joined['hour'].isin([10,11,12,13])]
afternoon = joined[joined['hour'].isin([14,15,16,17])]
evening = joined[joined['hour'].isin([18,19,20,21])]
night = joined[joined['hour'].isin([22,23,24,0,1,2,3,4,5])]

# print message
print('[INFO] - Combined data time of day counts:\n')
print('Morning sentences: ' + str(len(morning)))
print('Morning posts: ' + str(len(morning['id'].value_counts())))
print('Morning users: ' + str(len(morning['userid'].value_counts())))

print('\n')
print('Noon sentences: ' + str(len(noon)))
print('Noon posts: ' + str(len(noon['id'].value_counts())))
print('Noon users: ' + str(len(noon['userid'].value_counts())))

print('\n')
print('Afternoon sentences: ' + str(len(afternoon)))
print('Afternoon posts: ' + str(len(afternoon['id'].value_counts())))
print('Afternoon users: ' + str(len(afternoon['userid'].value_counts())))

print('\n')
print('Evening sentences: ' + str(len(evening)))
print('Evening posts: ' + str(len(evening['id'].value_counts())))
print('Evening users: ' + str(len(evening['userid'].value_counts())))

print('\n')
print('Night sentences: ' + str(len(night)))
print('Night posts: ' + str(len(night['id'].value_counts())))
print('Night users: ' + str(len(night['userid'].value_counts())))

# save to dataframe
morning.to_pickle('/home/waeiski/GIS/maphel_thirdplace/some_combined/combined/comb2015_morning.pkl')
noon.to_pickle('/home/waeiski/GIS/maphel_thirdplace/some_combined/combined/comb2015_noon.pkl')
afternoon.to_pickle('/home/waeiski/GIS/maphel_thirdplace/some_combined/combined/comb2015_afternoon.pkl')
evening.to_pickle('/home/waeiski/GIS/maphel_thirdplace/some_combined/combined/comb2015_evening.pkl')
night.to_pickle('/home/waeiski/GIS/maphel_thirdplace/some_combined/combined/comb2015_night.pkl')

# save to geopackage
joined.to_file('/home/waeiski/GIS/maphel_thirdplace/some_combined/combined/twinsta_2015_langid.gpkg', driver='GPKG')
