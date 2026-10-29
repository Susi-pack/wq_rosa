Water quality simulations for the paper led by Rosa.

[See paper here](https://docs.google.com/document/d/1jWefnUYo681QVpAxldnVLJbVBvlaG7_vRVIawxhQDVQ/edit?usp=sharing)

## How to use this repo
This is a project under `susi_26/projects/`.
Clone this repo under that folder and you should be able to run it!

(Don't worry: there should be no conflicts with nesting git repos.)

## TODO
- Prepare simulation scenarios
- Set the ash fertilization parameters
- Fix ditch shallowing. [See GitHub issue for making this parameter easier for users](https://github.com/Susi-pack/susi_26/issues/361)
- Double check the `drain_age` parameter, which is used to commpute the thick of the hummus layer, `h_mor`, but, is quite surprisingly, not used in the ditch shallowing!

## Simulation choices
This page is intended as a shared reference for the simulations.
Every choice made in the simulation should be reflected here.

- Metsäkeskus data filtering
    - When present, both dominant and subdominant canopylayers are simulated
    - There are stands with spruce in sfc 4 and 5. They are removed from the simulation
    - Only pine and spruce are allowed to be the dominant species
    - Removes spruce stands in fertility class {2,3}; pines can be in any sfc
- DNM
    - Not using DNM = -0.90; using DNM = -0.85 instead.
    - The default DNM is -0.60 cm

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
Still using a sample one.
Update this section when the true data comes.


