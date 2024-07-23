[![DOI](https://zenodo.org/badge/770497573.svg)](https://zenodo.org/doi/10.5281/zenodo.10810324)
# Biases in active land water storage capacity and its limitation of evapotranspiration in CMIP6 models
For full details on the methodology and results please refer to the manuscript folder above. For a step-by-step guide on how to reproduce the analysis, refer to the 'Instructions' below. All code is licensed under AGPL-v3, and the manuscript and data are licensed as CC-BY. Please review the individual directories and their LICENSE file for more information.

# Dry biases in land water storage and excessive soil moisture limitation in CMIP6 models
This repository provides code and intermediary data to reproduce the analysis of the abovementioned project. For full details on the methodology and results please refer to the 'manuscript' folder. For a step-by-step guide on how to reproduce the analysis, refer to the 'Instructions' below. 

All code is licensed under AGPL-v3, and the manuscript and data are licensed as CC-BY. Please review the individual directories and their LICENSE file for more information. You can cite the code in this repository as follows:

> Giardina et al. (2024). Biases in active land water storage capacity and its limitation of evapotranspiration in CMIP6 models: code and intermediary data. [https://zenodo.org/doi/10.5281/zenodo.10810324](https://zenodo.org/doi/10.5281/zenodo.10810324)


## Abstract
Accurate representation of plant water availability is crucial for climate modeling, due to its significant role in land-atmosphere interactions. Our study focuses on water storage dynamics and analyzes how soil moisture limitation is represented in Earth System Model (ESM) simulations of the Coupled Model Intercomparison Project phase 6 (CMIP6). We first quantify the long-term maximum annual depletion in water storage, contrasting model results with estimates based on satellite observations of terrestrial water storage from the Gravity Recovery and Climate Experiment (GRACE), as well as remotely sensed estimates of the water balance. Our analysis shows that CMIP6 models mostly underestimate the maximum annual soil moisture depletion, especially in the Amazon region. We further assess the frequency of soil moisture limitation in CMIP6 simulations against observations from solar-induced fluorescence (SIF) and GRACE, finding that ESMs generally overestimate this frequency. We evaluate our findings with analyses at 128 FLUXNET2015 sites, finding overall consistent results. Our study highlights the importance of improving the representation of plant water availability and land-atmosphere interactions in Earth System Models. Implementation of new model features could have large implications in predicting future climate on land.


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

