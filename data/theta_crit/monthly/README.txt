

Results produced using same years for GRACE and models (2007-2014 included), see Methods

files ending in _allscenarios: contains data from NINE CMIP6 models from both land-hist and historical scenarios (took all models available for all the variables used in this study)

OLD:
df_count_mrso.rds --> calculated with only 7 models in 2023!!! first version of the paper 


*** current df_count_GRACE and df_count_mrso_allscenarios use averaging to decrease the resolution of a dataset to match CMIP6, and nearest neighbour if it's a categorical variable. Corrected from before where bilinear interpolation was used for everything (default method when using resampling)



