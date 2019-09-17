## We rasterized the maps of population with access to improved water and improved sanitation from Pullan et al. (2014). The 1 km x 1km water.r and san.r raster objects should be first saved in data/ as "africa_wash_water_raster1km.rda" and "africa_wash_san_raster1km.rda," respectively. Additional data is available upon request.
## See original source at https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1001626 

## This code generates lists of targeted districts for water- and sanitation-based vaccine deployment strategies.

source("source/utils_mods.R")
reload_source()
options(error=recover)

## TOGGLE ME: baseline, high_campaignFreq, high_coverage, high_indirect, high_lifeExpect, high_vaccSupply, high_ve, low_campaignFreq, low_coverage, low_indirect, low_lifeExpect, low_vaccSupply, low_ve
sensitivityMod <- "baseline"

##############################
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
alloc_strategy_ls <- c("optimumWash", "optimumWat", "optimumSan") 
scenario_ls <- paste0(alloc_strategy_ls, "_who030718_v2")

method_coverage <- ifelse(sensitivityMod %in% c("baseline", "high_campaignFreq", "low_campaignFreq", "high_ve", "low_ve", "low_vaccSupply", "high_vaccSupply"), "wMedian1090 middle", ifelse(sensitivityMod == "high_coverage", "wMedian1090 ci_upper", ifelse(sensitivityMod == "low_coverage", "wMedian1090 ci_lower", NA))) 
rotate_every_nYears <- ifelse(sensitivityMod %in% c("baseline", "high_coverage", "low_coverage", "high_ve", "low_ve", "low_vaccSupply", "high_vaccSupply"), 3, ifelse(sensitivityMod == "high_campaignFreq", 2, ifelse(sensitivityMod == "low_campaignFreq", 5, NA))) 

global_input_fname <- paste0("source_inputs/vaccine_supply_", ifelse(sensitivityMod == "low_vaccSupply", "constant.csv", ifelse(sensitivityMod == "high_vaccSupply", "high.csv", "baseline.csv"))) 
all_inputs <- read_csv(global_input_fname) %>% dplyr::arrange(year)

countryLs <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GAB", "GHA", "GIN", "GMB", "GNB", "GBQ", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SDN", "SEN", "SLE", "SOM", "SSD", "SWZ", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE")

#######################
## Run each scenario 
fns <- substring(list.files(path = rasters_dir, pattern = "uWASH"), 1, 3)
countryLs2 <- countryLs[!(countryLs %in% fns)]
if(length(countryLs2)>0){
  wash_districts <- save_country_wash_raster(countryLs2, out_wd = rasters_dir) 
}


if(length(countryLs)==40){
  if(!file.exists(paste0(rasters_dir, "unimprovedWASH_2012_per_district.csv"))){
    fns <- list.files(path = rasters_dir, pattern = "uWASH")  
    wash_districts <- map_dfr(fns, function(fn){
      read_csv(paste0(rasters_dir, fn))
    })
    write_csv(wash_districts, paste0(rasters_dir, "unimprovedWASH_2012_per_district.csv"))
  } else{
    wash_districts <- read_csv(paste0(rasters_dir, "unimprovedWASH_2012_per_district.csv"))
  }
  

  for(k = 1:length(alloc_strategy_ls)){

    scen <- scenario_ls[k]
    cat(paste("Starting", scen), stdout())
    my_alloc_strategy <- alloc_strategy_ls[k]
   
    scen_targets <- identify_optimumWashPullan_targets(wash_districts, cntryCode_ls = countryLs, scenario = scen, vacc_inputs = all_inputs, method_coverage_prop = method_coverage, rotation = rotate_every_nYears, gen_in_wd = rasters_dir, gen_out_wd = targets_dir)  

    cat(paste("Ending",scen),stdout())

  }

}





