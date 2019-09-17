# The impact of geographic targeting of oral cholera vaccination in sub-Saharan Africa: a modeling study

This repository provides the source code and supplementary material for a modeling study on the impact of geographic targeting of oral cholera vaccination in sub-Saharan Africa. 

## Source

The source code for the primary models presented in the manuscript may be found in `source/`. The codes here demonstrate the specific settings and parameters in our models. Please note that the input data to run the source code are not posted in this repository, but descriptions for accessing this data exist at the top of each R script. 

Descriptions of the files are as follows:
  1. `generate_all_rasters_parallel.R`: Projects population at 5km x 5km scale for all study countries from 2018-2030. Crops and disaggregates incidence rate rasters for all study countries from Lessler and Moore et al., Lancet, 2018.
  1. `generate_wash_targets.R`: Processes 1km x 1km rasters on access to improved water and sanitation from Pullan et al., PLOS Medicine, 2014 into lists of targeted districts for water, sanitation, and WASH based vaccine deployment strategies for a given parameter set.
  1. `run_model_nonwash.R`: Run model for untargeted and burden based vaccine deployment strategies, for a given parameter set.
  1. `run_model_wash.R`: Run model for water and sanitation related vaccine deployment strategies, for a given parameter set.
  1. `process_all_models.R`: Processes non-wash and wash model outputs into non-cumulative and cumulative results (over study period years) and calculates confidence intervals for all parameter sets together.
  1. `make_main_plots.R`: Plot subpanels of main manuscript Figures 1 and 2.
  1. `make_main_plots_fig3.R`: Plot subpanels of main manuscript Figure 3.
  1. `utils_mods.R`: Utilities for targeting vaccine, modeling the health impacts of vaccination campaigns, and processing model output data. 
  1. `utils_figs.R`: Utilities for making the main manuscript figures.
  1. `utils_spatial.R`: Utilities for spatial- and raster-related functions employed by the model. Several of these functions will be replaced by a soon-to-be-released CRAN package in the future, as currently, there are missing file dependencies (last updated September 2019).

## Data

The input data files from external data sources may be found in `data/`. These files are called by the model codes or the supplement file `SM_geotargeting_ocv_ssa.Rmd`.

Descriptions of the files are as follows:
  * `WorldBank_GDP/`: GDP by country according to the WorldBank.
  * `ocv_ve_overtime`: Summary of vaccine efficacy data from Bi et al., Lancet Infectious Disease (2017).
  * `Review_OCV_campaign_coverage`: Literature review of oral cholera vaccination campaign coverage surveys.
  * `vacc_costs`: Literature review of oral cholera vaccination campaign cost surveys.
  * `who_cfrs`: Summary of case fatality ratios by country according to WHO Annual Cholera Report data from 1970-2016.
  * `whoannual_wppPop`: WHO Annual Cholera Report data by country and year from 1970-2017.
  * `WPP2017_SA4_MORT_F07_1_LIFE_EXPECTANCY_0_BOTH_SEXES`: Life expectancy data by country and year according to the United Nations World Population Prospects 2017 Revision.

## Source Inputs

We generated input data files for vaccine supply, and they may be found in `source_inputs/`.

## Country Targets

The `country_targets/` folder contains all of the district-level vaccine targets for each optimized vaccine deployment strategy according to the baseline parameter set in all study years. The files are organized by vaccine deployment strategy as follows:
  * `hrDistrict2`: case-logistics optimized
  * `hrDistrict2Incid`: rate-logistics optimized
  * `optimumDistrict`: case optimized
  * `optimumDistrictIncid`: rate optimized
  * `optimumSan`: sanitation optimized
  * `optimumWash`: WASH optimized
  * `optimumWat`: water optimized

## Generated Outputs

This folder contains the processed model outputs (from `process_all_models.R`) necessary for generating tables and figures in the supplement file `SM_geotargeting_ocv_ssa.Rmd`. 

The folders are named according to the sensitivity parameter set as follows:
  1. `baseline/`: Baseline parameter set reported in the main text.
  1. `high_campaignFreq/`: Vaccination campaigns every 2 years.
  1. `low_campaignFreq/`: Vaccination campaigns every 5 years.
  1. `high_coverage/`: 84\% vaccination campaign coverage. 
  1. `low_coverage/`: 50\% vaccination campaign coverage.
  1. `high_indirect/`: Upper bound vaccine indirect protection assumptions from Longini Jr. et al., PLOS Medicine 2007. 
  1. `low_indirect/`: No vaccine indirect protection.
  1. `high_lifeExpect/`: Population turnover rate of the inverse of 70 years (represents the lower bound of population turnover rate).
  1. `low_lifeExpect/`: Population turnover rate of the inverse of 56 yeras (represents the upper bound of population turnover rate).
  1. `high_vaccSupply/`: Vaccine supply increases linearly to a 95 million dose supply in 2030.
  1. `low_vaccSupply/`: Vaccine supply remains constant after 2019, achieving ad26 million dose supply in 2030.
  1. `high_ve/`: High vaccine efficacy assumptions.
  1. `low_ve/`: Low vaccine efficacy assumptions.

## Figures

This folder contains flowcharts that depict the rate optimized (`optimal_district_rate_flowchart`) and rate-logistics optimized (`hr_district_rate_flowchart`) vaccine deployment strategies. These images are pulled automatically into the supplement file `SM_geotargeting_ocv_ssa.Rmd`.



