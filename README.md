[![DOI](https://zenodo.org/badge/770497573.svg)](https://zenodo.org/doi/10.5281/zenodo.10810324)
# Large biases in the frequency of water limitation across Earth system models
This repository provides code and intermediary data to reproduce the analysis of the abovementioned project. For full details on the methodology and results please refer to the 'manuscript' folder. For a step-by-step guide on how to reproduce the analysis, refer to the 'Instructions' below.

All code is licensed under AGPL-v3, and the manuscript and data are licensed as CC-BY. Please review the individual directories and their LICENSE file for more information. You can cite the code in this repository as follows:
> Giardina et al. (2024). Biases in active land water storage capacity and its limitation of evapotranspiration in CMIP6 models: code and intermediary data. [https://zenodo.org/doi/10.5281/zenodo.10810324](https://zenodo.org/doi/10.5281/zenodo.10810324)


## Abstract
Water availability limits evapotranspiration on land, shaping the energy balance, land carbon uptake, and climate extremes. Despite its importance, Earth System Models struggle to capture where and how often water-limited conditions occur. Here we investigate the representation of water limitation and its link to land water storage capacity in simulations from the Sixth Coupled Model Intercomparison Project (CMIP6) driven by consistent observational atmospheric forcing. Using remotely sensed solar-induced vegetation fluorescence and terrestrial water storage, together with ecosystem flux observations, we find that CMIP6 models overestimate the frequency of water limitation by 14% over land and 26% in the tropics. Model overestimation occurs over 58% of the land area, and 78% in the tropics. These too frequent water-limited conditions are not conclusively linked to a potential underestimation of land water storage capacity in the models, hinting at gaps in how ESMs represent rooting depths, plant water uptake, and plant water-use strategies. Our study highlights the need for model development in these areas, with implications for projections of future climate on land.


## Instructions
First, clone this repo to your local computer:

```
git clone https://github.com/fgiardin/water_biases_CMIP6
```

To reproduce the analysis and figures, you can follow the steps described in the `analysis` folder. To avoid overwriting the dataframes already loaded in this repo, all scripts will generally save the output in the main directory of the project (aka the directory where this README also is). 

Below is an overview of the content of this repo. For detailed instructions, please consult the README files located in each subdirectory, or refer to the opening line of each script for a description of its purpose.

* `R`: contains all the R functions used in the analysis.
* `data-raw`: contains raw data and the scripts used to download, extract and process raw data. Very big data are not uploaded to this repo; please refer to the "Data availability" section of the manuscript to download them.
* `data`: contains processed data and output of the analysis.
* `manuscript`: contains manuscript and figures.

