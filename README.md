Water quality simulations for the paper led by Rosa.

[See paper here](https://docs.google.com/document/d/1jWefnUYo681QVpAxldnVLJbVBvlaG7_vRVIawxhQDVQ/edit?usp=sharing)

## How to use this repo
This is a project under `susi_26/projects/`.
Clone this repo under that folder and you should be able to run it!

(Don't worry: there should be no conflicts with nesting git repos.)

## TODO
- Implement thinning scenarios
- Improve coding style of run.py
- Double check ditch shallowing. [See GitHub issue for making this parameter easier for users](https://github.com/Susi-pack/susi_26/issues/361)
- Double check the `drain_age` parameter, which is used to commpute the thick of the hummus layer, `h_mor`, but, is quite surprisingly, not used in the ditch shallowing!
- Fix photo_parameters: change for each area.

## Simulation choices
This page is intended as a shared reference for the simulations.
Every choice made in the simulation should be reflected here.

- Metsäkeskus data filtering
    - When present, both dominant and subdominant canopylayers are simulated
    - There are stands with spruce in sfc 4 and 5. They are removed from the simulation
    - Only pine and spruce are allowed to be the dominant species
    - Removes spruce stands in fertility class {4,5}; pines can be in any sfc
- Scenarios:
    - See the `create_scenarios()` function for all the ditch and fertilization choices.

## How to get the missing data
Some files are too large to be tracked by GitHub, and they are gitignored.
Here's a list:
- Metsäkeskus .gpkg data used to create the allometry files.
- SUSI output results

### Metsäkeskus data
(NOTE: this is only useful if you wish to re-generate the allometry files.
If you only want to run SUSI, you don't need these files, because all the information is already contained in the allometry files).
Go to [https://avoin.metsakeskus.fi/aineistot/MV/Maakunta/](https://avoin.metsakeskus.fi/aineistot/MV/Maakunta/) and download the datasets for Pirkanmaa and Pohjois Pohjanmaa.

### Weather data
Both files hold daily weather for 1975-2025, in the same format as `CFw.csv`.
They were built from the FMI 10 km daily grids listed on [Paituli](https://paituli.csc.fi/).
The grids can be downloaded from the [Funet mirror](https://www.nic.funet.fi/index/geodata/ilmatiede/), as one NetCDF file per year and variable.
The variables used are daily mean, maximum and minimum temperature, precipitation, global radiation and relative humidity.

Each file holds a single grid cell, the one nearest the mean centre of the region's stands:
- Pirkanmaa: easting 318690, northing 6861775 (Ikaalinen)
- Pohjois-Pohjanmaa: easting 448690, northing 7171775 (Tyrnävä)

Coordinates are ETRS-TM35FIN metres, and go in the `longitude` and `latitude` columns.

Two values are not read directly from the grids:
- `hpa`: This is important, and points to an issue with the current `read_FMI_weather()`: in the FMI Paituli data there is no vapour pressure grid, but `read_FMI_weather()` needs it. The funny thing is that `read_FMI_weather()` actually requires that in order to derive RH, but RH is already in the data! So here we used the inverse formula to the one in the `read_FMI_weather()`, so that both cancel out. But this should be fixed properly later. This was recorded in [GH issue #365](https://github.com/Susi-pack/susi_26/issues/365). 
- `radiation` on 2016-12-30 and 2020-11-02: these days are missing from the radiation grids, so they are linearly interpolated from the neighbouring days.


