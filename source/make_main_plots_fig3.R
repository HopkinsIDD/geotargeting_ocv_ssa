## must run process_all_models.R first
## This code makes figure 3 from the main manuscript, which is based on all sensitivity parameter model runs.

source("source/utils_ms.R")
reload_source()
source("manuscripts/GAVI Impact Estimation/source/utils_figs.R")
options(error=recover)

alloc_strategy_ls <- c("optimumDistrict", "allDistrict", "hrDistrict2", "hrDistrict2Incid", "optimumDistrictIncid", "optimumWat", "optimumSan", "optimumWash") 
scenario_ls <- paste0(alloc_strategy_ls, "_who030718_v2")

fig_wd <- "figures"
gen_out_wd <- "generated_outputs"
sensitivity_ls <- c("baseline", "low_ve", "high_ve", "low_campaignFreq", "high_campaignFreq", "low_coverage", "high_coverage", "low_vaccSupply", "high_vaccSupply", "low_lifeExpect", "high_lifeExpect", "low_indirect", "high_indirect") 

dir.create(fig_wd, showWarnings = FALSE)

## import and process summary data
cum_summ <- get_cum_sensitivityRun_summaries(gen_output_dir = gen_out_wd, sensLs = sensitivity_ls)

tornadoDat <- plot_tornado_costper_oneway(cum_ssa = cum_summ, scenarioLs = scenario_ls, sensLs = sensitivity_ls, meas = "dalysAverted", alloc_strategy_bl = "hrDistrict2Incid", fig_out_wd = fig_wd)
