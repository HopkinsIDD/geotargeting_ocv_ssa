## must run generate_all_rasters_parallel first

source("source/utils_mods.R")
reload_source()

library(foreach)
library(doParallel)
cl <- makeCluster(5, outfile="target_model_nonwash_log.txt")
registerDoParallel(cl)
options(error=recover)

## TOGGLE ME: baseline, high_campaignFreq, high_coverage, high_indirect, high_lifeExpect, high_vaccSupply, high_ve, low_campaignFreq, low_coverage, low_indirect, low_lifeExpect, low_vaccSupply, low_ve
sensitivityMod <- "baseline" 

##############################################
## Set paths
rasters_dir <- "country_outputs/"
output_dir <- paste0("generated_outputs/", sensitivityMod, "/")
dir.create(output_dir, showWarnings = FALSE)
cf_dir <- "generated_outputs/cf/"
dir.create(cf_dir, showWarnings = FALSE)
targets_dir <- paste0("country_targets/", ifelse(sensitivityMod %in% c("baseline", "high_campaignFreq", "low_campaignFreq", "high_coverage", "low_coverage", "high_vaccSupply", "low_vaccSupply"), sensitivityMod, "baseline"), "/")
dir.create(targets_dir, showWarnings = FALSE)

###############################
## Set inputs for model runs 
alloc_strategy_ls <- c("allDistrict", "hrDistrict2Incid", "optimumDistrictIncid", "hrDistrict2", "optimumDistrict") 
scenario_ls <- paste0(alloc_strategy_ls, "_who030718_v2")

## settings for coverage, campaign frequency, VE, indirect effects, life expectancy, vaccine supply models
method_coverage <- ifelse(sensitivityMod %in% c("baseline", "high_campaignFreq", "low_campaignFreq", "high_ve", "low_ve"), "wMedian1090 middle", ifelse(sensitivityMod == "high_coverage", "wMedian1090 ci_upper", ifelse(sensitivityMod == "low_coverage", "wMedian1090 ci_lower", NA))) 
rotation <- ifelse(sensitivityMod %in% c("baseline", "high_coverage", "low_coverage", "high_ve", "low_ve"), 3, ifelse(sensitivityMod == "high_campaignFreq", 2, ifelse(sensitivityMod == "low_campaignFreq", 5, NA))) 
ve_scenario <- ifelse(sensitivityMod %in% c("baseline", "high_coverage", "low_coverage", "high_campaignFreq", "low_campaignFreq"), "base", ifelse(sensitivityMod == "high_ve", "high", ifelse(sensitivityMod == "low_ve", "low", NA))) 
indirect_fxn <- ifelse(sensitivityMod == "high_indirect", generate_indirect_incidence_sensitivity_mult("upper"), ifelse(sensitivityMod == "low_indirect", generate_indirect_incidence_sensitivity_mult("lower"), generate_indirect_incidence_mult()))
lfexpect <- ifelse(sensitivityMod == "low_lifeExpect", 1/56, ifelse(sensitivityMod == "high_lifeExpect", 1/70, 1/65))

global_input_fname <- paste0("source_inputs/vaccine_supply_", ifelse(sensitivityMod == "low_vaccSupply", "constant.csv", ifelse(sensitivityMod == "high_vaccSupply", "high.csv", "baseline.csv"))) 
all_inputs <- read_csv(global_input_fname) %>% dplyr::arrange(year)

countryLs <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GAB", "GHA", "GIN", "GMB", "GNB", "GBQ", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SDN", "SEN", "SLE", "SOM", "SSD", "SWZ", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE")

#######################
## Run each scenario 
foreach (k = 1:length(alloc_strategy_ls), 
  .export=c("countryLs","alloc_strategy_ls","scenario_ls","global_input_fname","generate_pct_protect_function","all_inputs","generate_flatline_multiplier")) %dopar% {

  reload_source()

  scen <- scenario_ls[k]
  cat(paste("Starting",scen), stdout())
  my_alloc_strategy <- alloc_strategy_ls[k]

  if(sensitivityMod %in% c("high_campaignFreq", "high_coverage", "high_vaccSupply", "low_campaignFreq", "low_coverage", "low_vaccSupply")){ ## all other sensitivity analyses use the same vaccine targets as the baseline parameter set

    ## identifies vaccination targets for all strategies except those related to access to water and sanitation: apply to all countries at once
    targets <- run_targeting_strategies(cntryCode_ls = countryLs,
                                       alloc_strategy = my_alloc_strategy,
                                       project_global_vacc = FALSE,
                                       scenario = scen,
                                       method_coverage_prop = method_coverage,
                                       rotate_every_nYears = rotation,
                                       targets_out_wd = targets_dir,
                                       rasters_out_wd = rasters_dir,
                                       global_vacc_fname = global_input_fname)

  }
  

  for (i in 1:length(countryLs)){
    cntry <- countryLs[i]
    outfile <- run_country_scenario(country = cntry,
                                    scenario = scen,
                                    alloc_strategy = my_alloc_strategy,
				                            target_inc_threshold = 1/1000,
				                            mu = lfexpect,
                                    ve_direct = generate_pct_protect_function(my_trunc_year=5, my_ve_scen = ve_scenario),
                                    years = all_inputs$year,
				                            indirect_mult = indirect_fxn,
                                    secular_trend_multiplier = generate_flatline_multiplier(),
                        				    use_partial_cover = FALSE, track_vac = TRUE,
                        				    out_wd = rasters_dir,
                                    cf_wd = cf_dir,
                        				    targets_wd = targets_dir,
                        				    gen_out_wd = output_dir)
    
  }

  cat(paste("Ending",scen),stdout())
}







