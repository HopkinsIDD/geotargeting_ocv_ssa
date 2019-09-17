## must put 2000, 2005, 2010, 2015, and 2020 WorldPop AFR .tif files in data/ folder first 
## Download from https://www.worldpop.org/geodata/summary?id=139
## must put mean annual cholera incidence rate rasters (rate.tif) from Lessler and Moore 2018 in data/folder second 
## Download from http://www.iddynamics.jhsph.edu/projects/cholera-dynamics

## this script generates all the country rasters needed for the analyses

rm(list=ls())
options(error=recover)

source("source/utils_ms.R")
reload_source()
source("source/taxdat_utils.R")

##### model settings #####
rasters_dir <- "country_outputs/"
dir.create(rasters_dir, showWarnings = FALSE)

disaggCoef = 4
countryLs <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GAB", "GHA", "GIN", "GMB", "GNB", "GBQ", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SDN", "SEN", "SLE", "SOM", "SSD", "SWZ", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE")

study_years <- 2018:2030
foreach (co = countryLs, .export=c("countryLs","study_years","rasters_dir","get_lambdas_stack","disaggCoef")) %dopar% {
    
    reload_source()
    cat(sprintf("Generating %s rasters, please be patient.\n",co))

    if (lookup_WorldPop_region(co) == "AFR") {
        save_country_rasters(co, years = study_years,
                             orig_continent_lambdas_stack = get_lambdas_stack(loc="africa"),
                             out_wd = rasters_dir,
                             disagg = disaggCoef,
                             partial_cover = FALSE,
                             partial_cover_method = "raster")
    } else {
        print(paste("cannot identify the WorldPop region for", co))
    }
}
