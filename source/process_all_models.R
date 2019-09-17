## This code processes the model outputs (run_model_wash.R and run_model_nonwash.R) into cumulative results and calculates confidence intervals

source("source/utils_mods.R")
reload_source()
options(error=recover)

## TOGGLE: baseline, high_campaignFreq, high_coverage, high_indirect, high_lifeExpect, high_vaccSupply, high_ve, low_campaignFreq, low_coverage, low_indirect, low_lifeExpect, low_vaccSupply, low_ve
sensitivityMod <- "baseline" 

##############################
## list of scenarios
alloc_strategy_ls <- c("allDistrict", "hrDistrict2Incid", "optimumDistrictIncid", "hrDistrict2", "optimumDistrict", "optimumWat", "optimumWash", "optimumSan") 
scenario_ls <- paste0(alloc_strategy_ls, "_who030718_v2")
health_outcomes <- c("cases", "casesAverted", "deathsAverted", "dalysAverted") 

c_out_wd <- paste0("country_outputs/")
gen_out_wd <- paste0("generated_outputs/", sensitivityMod)
cf_genout_wd <- "generated_outputs/cf"
targets_wd <- paste0("country_targets/", ifelse(sensitivityMod %in% c("baseline", "high_campaignFreq", "low_campaignFreq", "high_coverage", "low_coverage", "low_vaccSupply", "high_vaccSupply"), sensitivityMod, "baseline"))

dir.create(fig_wd, showWarnings = FALSE)

##############################
## process outputs
runsamps <- get_run_samples(scenarioLs = scenario_ls, gen_output_dir = gen_out_wd, cf_wd = cf_genout_wd)
write_csv(runsamps, paste0(gen_out_wd, "/master_processed_samples.csv"))
cumrunsamps <- process_cum_run_samples(runsamps)
write_csv(cumrunsamps, paste0(gen_out_wd, "/master_processed_cum_samples.csv"))

## calculate confidence intervals for all countries
runCI <- process_CI_ssa(runsamps, scenario_ls)
write_csv(runCI, paste0(gen_out_wd, "/master_processed_ci_ssa.csv"))
cumrunCI <- process_CI_ssa(cumrunsamps, scenario_ls)
write_csv(cumrunCI, paste0(gen_out_wd, "/master_processed_cum_ci_ssa.csv"))

## calculate confidence intervals for cumulative percentage of cases averted
percCACI <- process_CI_percCA(cumrunsamps, scenario_ls)
write_csv(percCACI, paste0(gen_out_wd, "/master_processed_ci_percCA.csv"))

## calculate confidence intervals for cost per DALY averted
costperCI <- process_CI_costper(cumrunsamps, sensitivityMod) 
write_csv(costperCI, paste0(gen_out_wd, "/master_processed_ci_costper.csv"))

## identify vaccination campaign targets
cntrytargets <- get_target_summary(targets_dir = targets_wd)
write_csv(cntrytargets, paste0(gen_out_wd, "/master_countryTargets_summary.csv"))
dtargets <- get_target_district_summary(targets_dir = targets_wd)
write_csv(dtargets, paste0(gen_out_wd, "/master_districtTargets_summary.csv"))

