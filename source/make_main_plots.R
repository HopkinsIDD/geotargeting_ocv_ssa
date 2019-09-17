## must run process_all_models.R first
## This code makes figures 1 and 2 from the main manuscript for a single parameter set (plus additional exploratory figures displayed for the baseline parameter set in the supplement).


source("source/utils_mods.R")
reload_source()
source("source/utils_figs.R")
options(error=recover)

## TOGGLE ME: baseline, high_campaignFreq, high_coverage, high_indirect, high_lifeExpect, high_vaccSupply, high_ve, low_campaignFreq, low_coverage, low_indirect, low_lifeExpect, low_vaccSupply, low_ve
sensitivityMod <- "baseline" 

## list of scenarios
alloc_strategy_ls <- c("allDistrict", "hrDistrict2Incid", "optimumDistrictIncid", "hrDistrict2", "optimumDistrict", "optimumWat", "optimumWash", "optimumSan") ## in order, the codes represent the following strategies: no targeting, rate-logistics optimized, rate optimized, case-logistics optimized, case optimized, water optimized, watsan optimized, sanitation optimized
scenario_ls <- paste0(alloc_strategy_ls, "_who030718_v2")
health_outcomes <- c("cases", "casesAverted", "deathsAverted", "dalysAverted") 

c_out_wd <- paste0("country_outputs/")
gen_out_wd <- paste0("generated_outputs/", sensitivityMod)
cf_genout_wd <- "generated_outputs/cf"
targets_wd <- paste0("country_targets/", ifelse(sensitivityMod %in% c("baseline", "high_campaignFreq", "low_campaignFreq", "high_coverage", "low_coverage", "low_vaccSupply", "high_vaccSupply"), sensitivityMod, "baseline"))
fig_wd <- paste0("figures/", sensitivityMod) 
dir.create(fig_wd, showWarnings = FALSE)


## read processed outputs
runsamps <- read_csv(paste0(gen_out_wd, "/master_processed_samples.csv"))
cumrunsamps <- read_csv(paste0(gen_out_wd, "/master_processed_cum_samples.csv"))
runCI <- read_csv(paste0(gen_out_wd, "/master_processed_ci_ssa.csv"))
cumrunCI <- read_csv(paste0(gen_out_wd, "/master_processed_cum_ci_ssa.csv"))
cntrytargets <- read_csv(paste0(gen_out_wd, "/master_countryTargets_summary.csv"))
dtargets <- read_csv(paste0(gen_out_wd, "/master_districtTargets_summary.csv"))
fvps <- cntrytargets %>% dplyr::filter(scenario == "hrDistrict2_who030718_v2") %>% group_by(year) %>% summarise(fvp = sum(fvp))


#### MAKE PLOTS ####

## Country targeting example (DR Congo)
target_cntry <- "COD"

## Figure 1A
cntry_lambda <- plot_choro_oneCountry_lambda(cntry = target_cntry, out_wd = c_out_wd, fig_out_wd = fig_wd, year = 2018)

## Figures 1B and 1C
for (strat in alloc_strategy_ls){
  do.call(plot_choro_oneCountry_targets,
          list(cntry = target_cntry, alloc_strategy = strat, fig_out_wd = fig_wd, years = 2018:2021, dtargets_summ = dtargets))
}

## Figures 2A and 2B
for (strat in alloc_strategy_ls){
  do.call(plot_choro_admin2Africa_fvp,
          list(dtargets_df = dtargets, allocation_strategy = strat, fig_out_wd = fig_wd)) ## default = allyrs
}

for (hOut in health_outcomes){
  ## Figure 2C
  do.call(plot_ts_cum_healthOutcomes_byScenario,
          list(gathered_df = cumrunCI, healthOutcome = hOut, scenarioLs = scenario_ls, fig_out_wd = fig_wd))
}
 
## Figure 2C
do.call(plot_ts_relativeCases_byScenario,
          list(gathered_df = cumrunCI, scenarioLs = scenario_ls, fig_out_wd = fig_wd, cumulative = TRUE)) 

  

  

  


  