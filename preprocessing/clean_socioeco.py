#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 10:37:21 2021

This script cleans up the socio-economic data for OLS and SLM analysis.

The grid database RTK contains obfuscated entries for cells that have less than 10 people and
the obfuscated value is -1, which are replaced by None in this script.

@author: Tuomas Väisänen
"""

import statistics
import geopandas as gpd
from sklearn.linear_model import LinearRegression


# read socio-economic data in
df = gpd.read_file('rtk2015_250_hma.gpkg')

# calculate kids and pensioners
df['kids_u18'] = df['he_0_2'] + df['he_3_6'] + df['he_7_12'] + df['he_13_15'] + df['he_16_17']
df['pensioners'] = df['he_65_69'] + df['he_70_74'] + df['he_75_79'] + df['he_80_84'] + df['he_85_']


# replace na values with mean income
for i, row in df.iterrows():
    if row['hr_ktu'] == -1:
        # remove obfuscation
        df.at[i, 'mean_income'] = None
    else:
        df.at[i, 'mean_income'] = row['hr_ktu']

# clean up income group data
for i, row in df.iterrows():
    if row['hr_pi_tul'] == -1:
        # remove obfuscation
        df.at[i, 'low_income'] = None
        df.at[i, 'med_income'] = None
        df.at[i, 'hi_income'] = None
    else:
        df.at[i, 'low_income'] = row['hr_pi_tul']
        df.at[i, 'med_income'] = row['hr_ke_tul']
        df.at[i, 'hi_income'] = row['hr_hy_tul']

# calculate income group relational size
df['low_income_prop'] = df['low_income'] / df['hr_tuy']
df['med_income_prop'] = df['med_income'] / df['hr_tuy']
df['hi_income_prop'] = df['hi_income'] / df['hr_tuy']

# calculate employment stats
for i, row in df.iterrows():
    if ((row['pt_tyovy'] <= 9) and (row['pt_tyovy'] > 0)):
        # remove obfuscation
        df.at[i, 'employed'] = None
        df.at[i, 'unemployed'] = None
    else:
        df.at[i, 'employed'] = row['pt_tyoll']
        df.at[i, 'unemployed'] = row['pt_tyott']

# calculate proportional employedness
df['employed_prop'] = df['employed'] / df['pt_tyovy']
df['unemployed_prop'] = df['unemployed'] / df['pt_tyovy']

# get outside workforce
df['not_in_workforce'] = df['pt_tyovu']

# calculate students
for i, row in df.iterrows():
    if ((row['pt_tyovu'] <= 9) and (row['pt_tyovu'] > 0)):
        # remove obfuscation
        df.at[i, 'students'] = None
    else:
        df.at[i, 'students'] = row['pt_opisk']
        
# calculate proportional students
df['students_prop'] = df['students'] / df['pt_vakiy']

# clean up education level data
for i, row in df.iterrows():
    if ((row['ko_ika18y'] <= 9) and (row['ko_ika18y'] > 0)):
        # remove obfuscation
        df.at[i, 'basic_ed'] = None
        df.at[i, 'high_school'] = None
        df.at[i, 'vocational_ed'] = None
        df.at[i, 'bachelor_ed'] = None
        df.at[i, 'masters_ed'] = None
    else:
        df.at[i, 'basic_ed'] = row['ko_perus']
        df.at[i, 'high_school'] = row['ko_yliop']
        df.at[i, 'vocational_ed'] = row['ko_ammat']
        df.at[i, 'bachelor_ed'] = row['ko_al_kork']
        df.at[i, 'masters_ed'] = row['ko_yl_kork']

# calculate education group relational size
df['basic_ed_prop'] = df['basic_ed'] / df['ko_ika18y']
df['high_school_prop'] = df['high_school'] / df['ko_ika18y']
df['vocational_ed_prop'] = df['vocational_ed'] / df['ko_ika18y']
df['bachelor_ed_prop'] = df['bachelor_ed'] / df['ko_ika18y']
df['masters_ed_prop'] = df['masters_ed'] / df['ko_ika18y']
df['med_ed_prop'] = (df['high_school'] + df['vocational_ed']) / df['ko_ika18y']
df['high_ed_prop'] = (df['bachelor_ed'] + df['masters_ed']) / df['ko_ika18y']

# clarify job structure data
df['total_jobs'] = df['tp_tyopy']
df['education_jobs'] = df['tp_p_koul'] / df['tp_tyopy']
df['entertainment_jobs'] = df['tp_r_taid'] / df['tp_tyopy']
df['intorg_jobs'] = df['tp_u_kans'] / df['tp_tyopy']
df['science_jobs'] = df['tp_m_erik'] / df['tp_tyopy']
df['horeca_jobs'] = df['tp_i_majo'] / df['tp_tyopy']
df['retail_jobs'] = df['tp_g_kaup'] / df['tp_tyopy']
df['health_jobs'] = df['tp_q_terv'] / df['tp_tyopy']
df['infocom_jobs'] = df['tp_j_info'] / df['tp_tyopy']
df['edsci_jobs'] = (df['tp_m_erik'] + df['tp_p_koul']) / df['tp_tyopy']
df['retail_ent_jobs'] = (df['tp_g_kaup'] + df['tp_r_taid']) / df['tp_tyopy']

# clean up living space data
for i, row in df.iterrows():
    if row['te_as_valj'] == -1:
        # remove obfuscation
        df.at[i, 'mean_livingspace'] = None
    else:
        df.at[i, 'mean_livingspace'] = row['te_as_valj']
        
# clean up mean square meters per apartment
for i, row in df.iterrows():
    if row['ra_as_kpa'] == -1:
        # remove obfuscation
        df.at[i, 'avg_sq_m'] = None
    else:
        df.at[i, 'avg_sq_m'] = row['ra_as_kpa']
        
# clean owner/renter counts
for i, row in df.iterrows():
    if row['te_omis_as'] == -1:
        # remove obfuscation
        df.at[i, 'owners'] = None
    else:
        df.at[i, 'owners'] = row['te_omis_as']
for i, row in df.iterrows():
    if row['te_vuok_as'] == -1:
        # remove obfuscation
        df.at[i, 'renters'] = None
    else:
        df.at[i, 'renters'] = row['te_vuok_as']
        
# count owner/renter proportions
df['owner_household_prop'] = df['owners'] / df['te_taly']
df['rent_household_prop'] = df['renters'] / df['te_taly']

# calculate building type proportions
df['tot_buildings'] = df['ra_raky']
df['res_buildings'] = df['ra_asrak']
df['res_buildings_prop'] = df['ra_asrak'] / df['ra_raky']
df['oth_buildings'] = df['ra_muut']
df['oth_buildings_prop'] = df['ra_muut'] / df['ra_raky']

# count dwelling type proportions
df['small_dwel_prop'] = df['ra_pt_as'] / df['ra_asunn']
df['block_dwel_prop'] = df['ra_kt_as'] / df['ra_asunn']

# save select few
df = df[['KUNTA', 'euref_x', 'euref_y', 'id_nro', 'vuosi', 'he_vakiy', 'he_naiset',
         'he_miehet', 'he_kika','upper_income_prop', 'lower_income_prop', 'kids_u18',
         'pensioners', 'hr_ktu', 'hr_mtu', 'low_income', 'med_income', 'hi_income',
         'low_income_prop', 'med_income_prop', 'hi_income_prop', 'employed',
         'unemployed', 'employed_prop', 'unemployed_prop',
         'students', 'students_prop', 'basic_ed', 'high_school', 'vocational_ed',
         'bachelor_ed', 'masters_ed', 'basic_ed_prop', 'high_school_prop',
         'vocational_ed_prop', 'bachelor_ed_prop', 'masters_ed_prop', 'med_ed_prop',
         'high_ed_prop', 'total_jobs', 'retail_ent_jobs',
         'education_jobs', 'entertainment_jobs', 'intorg_jobs', 'science_jobs',
         'horeca_jobs', 'retail_jobs', 'edsci_jobs', 'health_jobs', 'infocom_jobs',
         'mean_livingspace', 'avg_sq_m', 'tot_buildings', 'res_buildings',
         'res_buildings_prop', 'oth_buildings', 'oth_buildings_prop', 'small_dwel_prop',
         'block_dwel_prop', 'owner_household_prop', 'rent_household_prop', 'geometry']]

# read in twinsta points
points = gpd.read_file('twinsta_2015_langid.gpkg')

# spatial join
joined = gpd.sjoin(df, points)

# group by and count post and users
counts = joined.groupby(['id_nro'])['userid', 'id'].nunique().reset_index().rename(columns={'userid':'usercount', 'id':'postcount'})
platforms = joined.groupby(['id_nro','platform'])['platform'].count().rename('count').reset_index()
platforms = platforms.pivot(index='id_nro', columns=['platform'], values='count').reset_index()

# join counts to socio-economic df
df = df.merge(counts, on='id_nro')

# save to geopackage
df.to_file('RTK_data_5.gpkg', driver='GPKG')
