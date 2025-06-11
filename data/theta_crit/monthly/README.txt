

Results produced using same years for GRACE and models (2007-2014 included), see Methods

files ending in _allscenarios: contains data from NINE CMIP6 models from both land-hist and historical scenarios (took all models available for all the variables used in this study)

OLD:
df_count_mrso.rds --> calculated with only 7 models in 2023!!! first version of the paper


*** current df_count_GRACE and df_count_mrso_allscenarios use averaging to decrease the resolution of a dataset to match CMIP6, and average + 50% filter for landcover.
Corrected from before where bilinear interpolation was used for everything (default method when using "resample" function from terra package) --> these results are still in "grace counts with bilinear interp (old)"

cmip6_daily_theta_crit_count.rds --> not affected by abovementioned issue (focus on flux locations only, no resampling)


