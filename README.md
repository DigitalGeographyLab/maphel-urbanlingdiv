# maphel-urbanlingdiv
This repository contains the scripts used in the article "Mapping urban linguistic diversity with social media and population register data" published in CEUS (hopefully).

### Data requirements

* Twitter and Instagram data from the Helsinki Metropolitan Area from the year 2015.
  * Instagram data is legacy data we had collected in 2016 before the API was closed down. The Twitter data can be collected with [tweetsearcher](https://github.com/DigitalGeographyLab/tweetsearcher).
* Statistics Finland 250 m grid database from year 2015
  * Apply from [here](https://www.stat.fi/tup/ruututietokanta/index_en.html)
* Individual-level first language information aggregated into 250 m grid from 2015
  * Apply from [here](https://www.stat.fi/tup/mikroaineistot/index_en.html)
* Dynamic population data from Helsinki Metropolitan Area by Bergroth et al. (2022) from [here](https://zenodo.org/record/4726996#.Yfk1e9-xWF4)

### Pre-analysis steps

* The language detection of the social media was done with fastText using scripts from [Hiippala et al. (2020)](https://github.com/DigitalGeographyLab/maphel-finlang)
* Linguistic diversity of register data was calculated in the Statistics Finland secure environment FIONA with a similar script to [neighborhood_diversities.py](preprocessing/neighborhood_diversities.py).

### Suggested running order of scripts

| Step | Script | Description | Input | Output |
| ---- | :----- | :---------- | :---- | :----- |
| 1 | [combine_twitter_insta.py](preprocessing/combine_twitter_insta.py) | Combines Twitter and Instagram data | Instagram and Twitter point features geopackage | Twitter-Instagram combined point features geopackage |
| 2 | [neighborhood_diversities.py](preprocessing/neighborhood_diversities.py) | Calcualtes linguistic diversities across times of day | Output from step 1 and grid database | Grid database with diversity metrics |
| 3 | [neighborhood_diversities_no_timeofday.py](preprocessing/neighborhood_diversities_no_timeofday.py) | Calcualtes linguistic diversities for social media as a whole | Output from step 1 and grid database | Grid database with diversity metrics |
| 4 | [clean_socioeco.py](preprocessing/clean_socioeco.py) | Cleans socio-economic grid database | Raw RTK database file | Cleaned grid database |
| 5 | [join_dynpop.py](preprocessing/clean_socioeco.py) | Joins dynamic population to output from step 4 | Dynamic population data and output 4 | Grid database with dynamic population |
| 6 | [join_regdiv_to_socioeco.py](preprocessing/join_regdiv_to_socioeco.py) | Joins diversity metrics in registry with grid database from output 5 | Register data with linguistic diversity metrics and output from step 5 | Grid database |
| 7 | [user_langprofiles.py](preprocessing/user_langprofiles.py) | Calculates social media user linguistic profiles | Output from step 1 | Latex-formatted table |
| 8 | [moran_cluster.py](statistics/stability_socioeco.py) | Calculates clusters in register and social media grid data | Outputs from steps 6 and 3 | Geopackage with clusters |
| 9 | [stability_socioeco.py](preprocessing/stability_socioeco.py) | Classifies social media clusters based on temporal stability | Output from step 8 | Geopackage with stability classficiations |
| 10 | [extract_high_clusters.py](preprocessing/extract_high_clusters.py) | Extracts significant high linguistic diversity clusters | Output from step 9 | Geopackage with high diversity clusters |
| 11 | [kde_plot.py](visualization/kde_plot.py) | Plots linguistic diversity across times of day and the register data | Outputs from steps 6 and 3 | PNG file |
| 12 | [regression_ols_timeofday.py](statistics/regression_ols_timeofday.py) | Performs the OLS regression analysis |Outputs from steps 6 and 3 | Model files, VIF dataframes, error plots |
| 13 | SLM regression in [GeoDA](https://geodacenter.github.io/) | Run SLM regression |  Outputs from steps 6 and 3 | SLM model summaries |
| 14 | [OLS_summaries.py](statistics/OLS_summaries.py) | Prints OLS summaries | OLS model files from step 12 | Latex-formatted OLS summaries |
| 15 | [get_coef_table.py](statistics/get_coef_table.py) | Prints coefficient table from regression analyses | Output from step 14 | Latex-formatted table |
| 15 | [read_vif_files.py](statistics/read_vif_files.py) | Prints VIF test scroes | VIF dataframes from step 12 | Latex-formatted table |

### Note on analysis

The SLM analysis was conducted in [GeoDA](https://geodacenter.github.io/) with 1st order Queen contiguity neighborhoods.
