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

### Suggested running order of scripts

| Step | File | Description | Input | Output |
| ---- | :----- | :---------- | :---- | :----- |
| 0 | [place_cup.py](place_cup.py) | Places cup in correct position | Cup | Cup ready for coffee |
| 1 | [press_button.py](press_button.py) | Presses magic button | Cup ready for coffee | Cup full of coffee |
| 2 | [pick_cup.py](pick_cup.py) | Picks up cup full of coffee | Cup full of coffee | Coffee cup ready for drinking |
| 3 | [assess_flavor.r](assess_flavor.r) | Fits a linear model over flavor characteristics | Coffee cup ready for drinking | Happy researcher |
