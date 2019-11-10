
## utility functions 
## helper function to reload all packages and source code

reload_source <- function(){
    options(warn=1,error=traceback)
    source("source/utils_mods.R")
    library(sf)
    library(knitr)
    library(dplyr)
    library(broom)
    library(ggplot2)
    library(raster)
    library(rgdal)
    library(readxl)
    library(taxdat)
    library(purrr)
    library(readr)
    library(ggthemes)
    library(rworldmap)
    library(maps)
    library(RColorBrewer)
    library(gridExtra)
    library(readODS)
    library(tidyr)
}


#' @name get_lambdas_stack
#' @title get_lambdas_stack
#' @description Import mean annual incidence (2010-2016) rasters from Lessler and Moore (2018)
#' @param loc
#' @return R6 pointer to lambda stack 
get_lambdas_stack <- function(loc="africa"){
    if(tolower(loc) == "africa" | loc == "AFR"){
        lambdas_stack = R6raster$new("data/rate.tif")
    } 

    raster::projection(lambdas_stack$raster) = "+proj=longlat + ellps=WGS84 +datum=WGS84 +no_defs"
    return(lambdas_stack)
}


#' @name pop_weighted_wash_districts
#' @title pop_weighted_wash_districts
#' @description Aggregates Pullan et al. (2012) WASH coverage rasters at country-level to district according to population-weighted mean
#' @param cntry_code
#' @param cntry_wash_r6
#' @param cntry_pop_r6
#' @param sublvl
#' @param varname
#' @param partial_cover
#' @param partial_cover_method
#' @return 
pop_weighted_wash_districts <- function(cntry_code, cntry_wash_r6, cntry_pop_r6, cntry_pop_district, sublvl = 2, varname, partial_cover = FALSE, partial_cover_method = "raster"){
  
  ## multiply wash coverage by population (used as a weight)
  tmp_num <- R6raster$new(raster::overlay(cntry_wash_r6$raster, cntry_pop_r6$raster, fun = function(r1,r2){r1*r2}))
  tmp_num_district <- apply_to_all_sublevels(cntry_code, tmp_num,
                          ISO_level = sublvl, 
                          trim_small = 0.05,
                          partial_cover = partial_cover, 
                          method = partial_cover_method,
                          fun = function(rast){if(!all(is.na(rast))){
                                rc <- sum(rast, na.rm=TRUE)
                              } else {
                                rc <- NA
                              }
                          })

  ## summarise by values because districts are 
  num <- data.frame(district = names(tmp_num_district), val = tmp_num_district, row.names = NULL, stringsAsFactors = FALSE) %>%
    group_by(district) %>%
    summarise(val = sum(val)) %>% ungroup
  pop <- data.frame(district = names(cntry_pop_district), pop = cntry_pop_district, row.names = NULL, stringsAsFactors = FALSE) %>%
    group_by(district) %>%
    summarise(pop = sum(pop)) %>% ungroup

  new_df <- full_join(num, pop, by = c("district")) %>%
    dplyr::mutate(var = val/pop) %>%
    dplyr::select(district, var) %>%
    dplyr::rename(!!varname := var) 

  return(new_df)
}


#' @name save_country_wash_raster
#' @title save_country_wash_raster
#' @description Calculate probability of access to unimproved water, sanitation, and WASH and write to file
#' @param country_ls
#' @param year
#' @param orig_continent_wash_wd
#' @param out_wd
#' @param partial_cover
#' @param partial_cover_method
#' @return 
save_country_wash_raster <- function(country_ls,
                                year = 2012, ## year of Pullan et al wash coverage estimates
                                orig_continent_wash_wd = "data/",
                                out_wd = "country_outputs/",
                                partial_cover = FALSE,
                                partial_cover_method = "raster"
                                ){

  ## outputs from shifting origins
  load(paste0(orig_continent_wash_wd, "africa_wash_water_raster1km.rda")) ## water.r
  load(paste0(orig_continent_wash_wd, "africa_wash_san_raster1km.rda")) ## san.r

  ## process country wash data
  washDat <- map_dfr(country_ls, function(country){

    cat(sprintf("saving %s districts \n", country))

    ## set up shapefiles and case and incidence data
    cntry_code <- fix_country_name(country)
    who_region <- lookup_WorldPop_region(cntry_code)
    sublevel <- 2 

    cntry_shp <- tryCatch(get_country_shapefile(cntry_code),
                          error = function(e) {
                              return(NA) })

    if(class(cntry_shp)[1]=="logical"){
        warning(sprintf("not saving any files for %s since the country name isn't recognized", country))
        return(NA)
    }

    ## load pop data
    ## use a snap out version of the country raster to use when loading population estimates
    
    save_layer <- R6raster$new(raster(paste0("country_outputs/", cntry_code, "_lambda_per_cell.tif"))[[1]])
    pop_per_cell <- load_population_estimates(region = who_region, save_layer, year)
    ext <- extent(pop_per_cell$raster)
    rm(save_layer)
    gc()
    
    ## load wash data
    orig_cntry_wat_rast <- R6raster$new(raster::crop(water.r, ext))
    orig_cntry_san_rast <- R6raster$new(raster::crop(san.r, ext))
    cntry_wat_rast <- R6raster$new(raster::aggregate(orig_cntry_wat_rast$raster, 5))
    cntry_san_rast <- R6raster$new(raster::aggregate(orig_cntry_san_rast$raster, 5))
    rm(orig_cntry_wat_rast, orig_cntry_san_rast)
    gc()

    ## calculate unimproved wash coverage
    cntry_uWat_rast <- R6raster$new(calc(cntry_wat_rast$raster, fun = function(r1){1-r1}))
    cntry_uSan_rast <- R6raster$new(calc(cntry_san_rast$raster, fun = function(r1){1-r1}))
    cntry_uWash_rast <- R6raster$new(overlay(cntry_uWat_rast$raster, cntry_uSan_rast$raster, fun = function(r1,r2){return(r1+r2-(r1*r2))})) ## pr(unimproved water or unimproved sanitation)
    rm(cntry_wat_rast, cntry_san_rast)
    gc()

    ## aggregate population to districts as denominator
    pop_per_district <- apply_to_all_sublevels(cntry_code, pop_per_cell,
                          ISO_level = sublevel, 
                          trim_small = 0.05,
                          partial_cover = partial_cover, 
                          method = partial_cover_method,
                          fun = function(rast){if(!all(is.na(rast))){
                                rc <- sum(rast, na.rm=TRUE)
                              } else {
                                rc <- NA
                              }
                          })


    printall <- function(r1){
      print(deparse(substitute(r1)))
      print(r1)
      print(extent(r1))
      print(origin(r1))
      print("**********************")
    }
    printall(pop_per_cell$raster)
    printall(cntry_uWat_rast$raster)
    printall(cntry_uSan_rast$raster)
    printall(cntry_uWash_rast$raster)
    
    ## unimproved probabilities
    wat_per_district <- pop_weighted_wash_districts(cntry_code, cntry_uWat_rast, pop_per_cell, pop_per_district, sublvl = sublevel, varname = "wat")
    san_per_district <- pop_weighted_wash_districts(cntry_code, cntry_uSan_rast, pop_per_cell, pop_per_district, sublvl = sublevel, varname = "san")
    wash_per_district <- pop_weighted_wash_districts(cntry_code, cntry_uWash_rast, pop_per_cell, pop_per_district, sublvl = sublevel, varname = "wash")
    rm(cntry_uWat_rast, cntry_uSan_rast, cntry_uWash_rast, pop_per_cell)
    gc()

    ## rm duplicate districts
    ppd <- data.frame(district = names(pop_per_district), pop = pop_per_district, row.names = NULL, stringsAsFactors = FALSE) %>%
      group_by(district) %>%
      summarise(pop = sum(pop)) %>% ungroup

    full_df <- full_join(wat_per_district, san_per_district, by = c("district")) %>%
      full_join(wash_per_district, by = c("district")) %>%
      full_join(ppd, by = c("district")) %>%
      rowwise %>%
      dplyr::mutate(cntry_code = fix_country_name(unlist(strsplit(district, split = "_"))[1])) %>%
      ungroup %>% 
      dplyr::arrange(district) %>%
      dplyr::mutate(index = 1:nrow(.))
	  print(paste("adding full_df for", cntry_code))

    write_csv(full_df, paste0(out_wd, cntry_code, "_uWASH_2012_per_district.csv")) ## delete files after finished
    return(full_df)
  })

  return(washDat)
}


#' @name save_country_rasters
#' @title save_country_rasters
#' @description Imports incidence rate and population data for a given country and saves raster and district summaries for use by other parts of the model
#' @param country
#' @param years
#' @param orig_continent_lambdas_stack
#' @param out_wd
#' @param disagg
#' @param partial_cover
#' @param partial_cover_method
#' @return 
save_country_rasters <- function(country,
                                 years = 2018:2035,
                                 orig_continent_lambdas_stack = get_lambdas_stack(loc="africa"),
                                 out_wd = "country_outputs/",
                                 disagg = 4,
                                 partial_cover = FALSE,
                                 partial_cover_method = "raster"
                                 ){
    
    cat(sprintf("saving %s rasters \n", country))

    ## set up shapefiles and case and incidence data
    cntry_code <- fix_country_name(country)
    who_region <- lookup_WorldPop_region(cntry_code)
    sublevel <- 2 

    cntry_shp <- tryCatch(get_country_shapefile(cntry_code),
                          error = function(e) {
                              return(NA) })

    if(class(cntry_shp)[1]=="logical"){
        warning(sprintf("not saving any files for %s since the country name isn't recognized", country))
        return(NA)
    }

    if(!dir.exists(out_wd)){dir.create(out_wd)}

    ## extract lambda
    orig_cntry_lambdas_stack <- crop_raster_to_shapefile(orig_continent_lambdas_stack, cntry_shp)
    # lambda_per_cell <- orig_cntry_lambdas_stack
    lambda_per_cell_snapout <- R6raster$new(raster::disaggregate(orig_cntry_lambdas_stack$raster, disagg))
    save_layer <- R6raster$new(lambda_per_cell_snapout$raster[[1]])
    
    rm(orig_cntry_lambdas_stack, orig_continent_lambdas_stack)
    gc()
  
    if(!file.exists(paste0(out_wd, cntry_code, "_lambda_per_cell.tif"))){

      lambda_per_cell <- extract_country_from_raster(lambda_per_cell_snapout, cntry_shp, partial_cover = partial_cover, method = partial_cover_method, trim_small = 0.05)
      ## save lambda file
      writeRaster(lambda_per_cell$raster, paste0(out_wd, cntry_code, "_lambda_per_cell.tif"))

    } else{
      ## read existing saved lambda per cell
      lambda_per_cell <- R6raster$new(paste0(out_wd, cntry_code, "_lambda_per_cell.tif"))
    }
     
    
    ## calculate and save cases and population by year
    sapply(years, 
      function(year){
          print(sprintf("Processing %s %s objects", cntry_code, year))

          ## non-partial cover version
          tmp_pop <- load_population_estimates(region = who_region, save_layer, year)
         ## partial cover version
         pop_per_cell <- extract_country_from_raster(tmp_pop, cntry_shp, partial_cover = partial_cover, method = partial_cover_method, trim_small = 0.05)

         ## calculate cases per cell from lambdas and population (partial cover already applied to population)
         cases_per_cell <- R6raster$new(lambda_per_cell$raster * pop_per_cell$raster)

         nlayers <- dim(cases_per_cell$raster)[3]
         ## summarize data for district
         cases_per_district <- apply_to_all_sublevels(cntry_code,
                                 cases_per_cell,
                                 ISO_level = sublevel, 
                                 trim_small = 0.05,
                                 partial_cover = FALSE, ## partial cover already applied to cases through population multiplication
                                 method = partial_cover_method,
                                 fun = function(x){if(!all(is.na(x))){
                                       rc <- sum(x, na.rm=TRUE)
                                     } else {
                                       rc <- NA
                                     }
                                 })

         ## save and gc
         print(sprintf("Saving %s %s: pop_per_cell, cases_per_district", cntry_code, year))
         writeRaster(pop_per_cell$raster, paste0(out_wd, cntry_code, "_", year, "_pop_per_cell.tif"))
         saveRDS(cases_per_district, paste0(out_wd, cntry_code, "_", year,  "_cases_per_district.rds"))

         rm(cases_per_cell, pop_per_cell)
         gc()

          ## extracts population per district but then puts into a matrix to match the MCMCs
          pop_per_district <- rep(apply_to_all_sublevels(cntry_code,
                                    tmp_pop,
                                    ISO_level = sublevel,
                                    partial_cover = partial_cover, # if FALSE, some district pops are 0
                                    trim_small = 0.05,
                                    method = partial_cover_method,
                                    fun=function(x){if(!all(is.na(x))){
                                        rc <- sum(x, na.rm=TRUE)
                                      } else {
                                        rc <- NA
                                      }
            }), nlayers) %>%
                   matrix(.,
                          ncol=ncol(cases_per_district),
                          byrow = TRUE,
                          dimnames=list(1:nlayers,names(.)[1:ncol(cases_per_district)]))

          print(sprintf("Saving %s %s: pop_per_district", cntry_code, year))
          saveRDS(pop_per_district, paste0(out_wd, cntry_code, "_", year, "_pop_per_district.rds"))

          rm(pop_per_district, cases_per_district, tmp_pop)
          gc()
    })

    rm(lambda_per_cell, save_layer)
    gc()

    cat(sprintf("%s rasters saved \n", country))
    return(NA)
}


#' @name run_targeting_strategies
#' @title run_targeting_strategies
#' @description Identify vaccination strategies for all non-WASH-based vaccination deployment strategies
#' @param scenario
#' @param alloc_strategy
#' @param cntryCode_ls
#' @param project_global_vacc
#' @param rotate_every_nYears
#' @param targets_out_wd
#' @param global_vacc_fname
#' @return
run_targeting_strategies <- function(scenario,
                                    alloc_strategy,
                                    cntryCode_ls,
                                    project_global_vacc = FALSE,
                                    method_coverage_prop = "wMedian1090 middle", # string with 2 words: "wMeanSamp/wMedianSamp middle/ci_lower/ci_upper", NA
                                    rotate_every_nYears = 3,  
                                    rasters_out_wd = "country_outputs/",
                                    targets_out_wd = "country_targets/",
                                    global_vacc_fname){

  ## import file with year, global vaccine availability, and target coverage proportion; project forward to 2030 if needed
  if(project_global_vacc == TRUE){
    inData <- read_csv(global_vacc_fname) %>% dplyr::arrange(year)
    ## linear model
    fit <- lm(vaccines ~ year, data = inData)
    pred <- predict(fit, newdata = data.frame(year = inData$year))
    pred_cl <- ifelse(pred < 0, 0, pred)
    
    all_inputs <- tbl_df(data.frame(year = inData$year, origData = inData$vaccines, vaccines = signif(pred_cl, 3), coverage = inData$coverage[1])) %>%
      dplyr::filter(year >= 2019)

  } else {
      all_inputs <- read_csv(global_vacc_fname) %>% dplyr::arrange(year) %>%
        dplyr::filter(!is.na(vaccines))
  }

  ## get vaccination coverage level for the scenario
  if(!is.na(method_coverage_prop)){
      covEstimate <- get_vacc_coverage() %>% dplyr::filter(method == method_coverage_prop) %>% dplyr::select(estimate) %>% unlist
      all_inputs <- all_inputs %>%
        dplyr::mutate(coverage = covEstimate)
    } else{
      next
    }

  ## Assign targeting function according to the vaccine deployment strategy
  if (alloc_strategy == "optimumDistrictIncid"){ ## rate optimized
    print(paste("targeting by", alloc_strategy))
    my_targeting_function = identify_optimumDistrictIncid_targets
  } else if (alloc_strategy == "optimumDistrict"){ ## case optimized
    print(paste("targeting by", alloc_strategy))
    my_targeting_function = identify_optimumDistrict_targets
  } else if (alloc_strategy == "allDistrict"){ ## untargeted
    print(paste("targeting by", alloc_strategy))
    my_targeting_function = identify_allDistrict_targets 
  } else if (alloc_strategy == "hrDistrict2"){ ## case-logistics optimized
    print(paste("targeting by", alloc_strategy))
    my_targeting_function = identify_hrDistrict2_targets
  } else if (alloc_strategy == "hrDistrict2Incid"){ ## rate-logistics optimized
    print(paste("targeting by", alloc_strategy))
    my_targeting_function = identify_hrDistrict2Incid_targets
  } 
    
  for (i in 1:nrow(all_inputs)){
    if(all_inputs$vaccines[i] > 0 & all_inputs$coverage[i] > 0){
      targets <- do.call(my_targeting_function, 
                        list(cntryCode_ls = cntryCode_ls, 
                          scenario = scenario, 
                          rotate_every_nYears = rotate_every_nYears, 
                          global_vacc_avail = as.numeric(all_inputs$vaccines[i]), 
                          year = all_inputs$year[i], 
                          target_coverage_prop = as.numeric(all_inputs$coverage[i]), 
                          gen_out_wd = targets_out_wd, 
                          out_wd = rasters_out_wd))
    } else{
      print(paste("No vaccine allocation or 0% vaccination coverage in ", all_inputs$year[i]))
      next
    }
  }
    
  return()
}


#' @name identify_optimumDistrictIncid_targets
#' @title identify_optimumDistrictIncid_targets
#' @description Allocate vaccines with the rate optimized vaccination deployment strategy
#' @param cntryCode_ls
#' @param scenario
#' @param global_vacc_avail Total number of vaccines available to allocate, globally 
#' @param target_coverage_prop Vaccination coverage 
#' @param doses_per_person Two-dose full course
#' @param rotate_every_nYears
#' @param year
#' @param out_wd
#' @param gen_out_wd
#' @return list with two dataframes, one for district-level (fullDistrict) targets and one for cell-level (partialDistrict) targets
identify_optimumDistrictIncid_targets <- function(cntryCode_ls,
                              scenario,
                              global_vacc_avail,
                              target_coverage_prop,
                              doses_per_person = 2,
                              rotate_every_nYears = 3,
                              partial_cover = FALSE,
                              partial_cover_method = "raster",
                              year,
                              out_wd = "country_outputs/",
                              gen_out_wd = "country_targets/"
                              ){

    cat(sprintf("allocate vaccines with optimum district targeting by incidence for %s: %s \n", year, paste(cntryCode_ls, collapse=", ")))
    
    inputs <- data.frame(country = cntryCode_ls, year = year) %>% 
        dplyr::mutate(country_year = paste(country, year, sep="--")) %>%
        dplyr::select(country_year) %>% unname %>% unlist
    all_district_ranks <- dplyr::bind_rows(mapply(FUN = rank_districtsIncid_by_country, country_year=inputs, scenario=scenario, rotate_every_nYears=rotate_every_nYears, in_wd = out_wd, out_wd = gen_out_wd, SIMPLIFY = FALSE)) %>%
        dplyr::filter(!skip)
    
    ## allocate vaccine to districts with highest cases
    all_district_ranks <- tbl_df(all_district_ranks[order(all_district_ranks$incidence, decreasing = TRUE),]) %>%
        dplyr::mutate(id = seq_along(pop),
                    vacc_pop_coverage = pop * as.numeric(target_coverage_prop) * doses_per_person,
                    vacc_pop_cumsum = cumsum(vacc_pop_coverage),
                    vacc_used = ifelse(vacc_pop_cumsum <= global_vacc_avail, vacc_pop_coverage, 0),
                    fvp = vacc_used/doses_per_person,
                    fullDistrict_target = ifelse(vacc_used > 0, TRUE, FALSE))

    fullDistrict_target_names <- all_district_ranks$district[which(all_district_ranks$fullDistrict_target)]
    cat(paste0("Fully-targeted districts: ", paste(fullDistrict_target_names, collapse = ", "),"\n"))


    ## calculate number of vaccines remaining for district receiving partial allocation
    partialDistrict_vacc_pop_size <- global_vacc_avail - sum(all_district_ranks$vacc_used, na.rm = TRUE)
    ## identify highest ranked district that did not receive vaccine as the district that will receive a partial allocation
    ## NOTE: adding the criteria that pop needs to be greater than 0 to be targeted. This can occur with small polygons
    partialDistrict_target_id <- all_district_ranks$id[min(which(all_district_ranks$vacc_used==0 & all_district_ranks$pop>0))]
    partialDistrict_target_name <- all_district_ranks$district[partialDistrict_target_id]
    print(partialDistrict_vacc_pop_size)
    print(partialDistrict_target_id)
    cat(sprintf("Partially-targeted district: %s (%.1f vaccinated) \n ",
                partialDistrict_target_name, partialDistrict_vacc_pop_size))

    ## import cell-level rasters for the partial district country
    sublevel <- 2
    cntry_code <- all_district_ranks$cntry_code[partialDistrict_target_id]
    partialDistrict_target_poly_id <- all_district_ranks$index[partialDistrict_target_id] # index is the cntry-specific polygon ID
    cntry_code <- fix_country_name(cntry_code)
    lambda_per_cell <- R6raster$new(paste0(out_wd,cntry_code,"_lambda_per_cell.tif"))
    pop_per_cell <- R6raster$new(paste0(out_wd, cntry_code,"_",year, "_pop_per_cell.tif"))
    cases_per_cell <- R6raster$new(lambda_per_cell$raster * pop_per_cell$raster)
    cntry_shp <- get_country_sublevels(cntry_code, sublevel)

    ## get rasters for partially targeted district
    partialDistrict_lambda_rast <- extract_country_from_raster(lambda_per_cell, cntry_shp$geometry[partialDistrict_target_poly_id], partial_cover = FALSE, method = "raster")
    partialDistrict_pop_rast <- extract_country_from_raster(pop_per_cell, cntry_shp$geometry[partialDistrict_target_poly_id], partial_cover = partial_cover, method = partial_cover_method, trim_small = 0.05)
    partialDistrict_cases_rast <- R6raster$new(partialDistrict_lambda_rast$raster * partialDistrict_pop_rast$raster)
    
    partialDistrict_cases_summ_rast <- aggregate_raster_xlayers(partialDistrict_cases_rast, median, na.rm = TRUE)

    ## allocate vaccine to the grid cells with greatest incidence in partialDistrict at a given coverage level until there is no more
    partialDistrict_lambda_summ_rast <- aggregate_raster_xlayers(partialDistrict_lambda_rast, median, na.rm = TRUE)  

    ## some cells are NA, these ranks need to be last
    cell_rank <- data.frame(
        cntry_code = cntry_code,
        district = all_district_ranks$district[partialDistrict_target_id],
        year = year,
        index = order(values(partialDistrict_lambda_summ_rast$raster), decreasing=T, na.last=TRUE),
        incidence = values(partialDistrict_lambda_summ_rast$raster)[order(values(partialDistrict_lambda_summ_rast$raster), decreasing=T, na.last=TRUE)],
        cases = values(partialDistrict_cases_summ_rast$raster)[order(values(partialDistrict_lambda_summ_rast$raster), decreasing=T, na.last=TRUE)],
        pop = values(partialDistrict_pop_rast$raster)[order(values(partialDistrict_lambda_summ_rast$raster), decreasing=T, na.last=TRUE)])

    ## get all but the last cell to target here
    cell_rank <- cell_rank %>%
    dplyr::mutate(vacc_pop_coverage = pop * target_coverage_prop * doses_per_person,
                  vacc_pop_cumsum = cumsum(vacc_pop_coverage),
                  vacc_used = ifelse(vacc_pop_cumsum <= partialDistrict_vacc_pop_size, vacc_pop_coverage, 0),
                  fullCell_target = ifelse(vacc_used > 0, TRUE, FALSE))

    partialDistrict_target_cell_ids <- cell_rank[which(cell_rank$fullCell_target),]$index

    ## allocate last vaccines to highest-ranked untargeted cell, partially
    partialCell_vacc_pop_size <- partialDistrict_vacc_pop_size - sum(cell_rank$vacc_used, na.rm = TRUE)

    cat(paste0(" Full district vaccination: ", sum(all_district_ranks$vacc_used, na.rm = TRUE),"\n"))
    cat(paste0(" Partial district vaccination: ", partialDistrict_vacc_pop_size,"\n"))
    cat(paste0(" Total vaccination: ", sum(all_district_ranks$vacc_used, na.rm = TRUE) +
                                      partialDistrict_vacc_pop_size,"\n"))


    partialCell_target_cell_id <-cell_rank$index[min(which(cell_rank$vacc_used == 0), na.rm=T)] 

    if(!is.na(partialCell_target_cell_id)){
      cell_rank[which(cell_rank$index==partialCell_target_cell_id),]$vacc_used <- partialCell_vacc_pop_size
      partialCell_target_prop <- (partialCell_vacc_pop_size/doses_per_person)/cell_rank[which(cell_rank$index==partialCell_target_cell_id),]$pop

      print(paste("partial cell target ID", partialCell_target_cell_id))
      print(paste("vacc pop size", partialCell_vacc_pop_size))
      print(paste("target prop", partialCell_target_prop))
    } else{
      print(paste(partialCell_vacc_pop_size, "partial cell vaccines not allocated."))
    }

    cell_rank <- cell_rank %>%
      dplyr::mutate(fvp = vacc_used/doses_per_person,
                    prop_vaccinated = ifelse(pop==0, 0, fvp/pop))

    rm(lambda_per_cell, cases_per_cell, pop_per_cell)
    gc()

    if(!dir.exists(gen_out_wd)){
        dir.create(gen_out_wd)
    }

    ## Need these outputs to create vac_template_rast by country in a separate function
    write_csv(all_district_ranks, paste0(gen_out_wd, scenario, "_", year, "_districtTargets.csv"))
    write_csv(cell_rank, paste0(gen_out_wd, scenario, "_", year, "_cellTargets.csv"))

    return(list(district_targets = all_district_ranks, cell_targets = cell_rank))
}


#' @name rank_districtsIncid_by_country
#' @title rank_districtsIncid_by_country
#' @description Rank districts by incidence within a single country
#' @param country_year
#' @param scenario
#' @param rotate_every_nYears
#' @param in_wd
#' @param out_wd
#' @return dataframe with ordered districts and incidence by country
rank_districtsIncid_by_country <- function(country_year, 
                                  scenario,
                                  rotate_every_nYears = 3,
                                  in_wd = "country_outputs/",
                                  out_wd = "country_targets/"
                                  ){

    country <- unlist(strsplit(country_year, split="--"))[1]
    year <- as.numeric(unlist(strsplit(country_year, split="--"))[2])
    cat(sprintf("ranking %s districts in %s \n", country, year))

    ## get list of fully targeted districts in previous x years
    skipYears <- year-(1:(rotate_every_nYears-1))
    skip_ls <- skip_district_targets(skipYears, scenario, out_wd)
    
    cntry <- fix_country_name(country)
    ## import district-level data
    cases_per_district_samps <- readRDS(file = paste0(in_wd,cntry,"_",year,"_cases_per_district.rds"))
    pop_per_district <- readRDS(file = paste0(in_wd,cntry,"_",year,"_pop_per_district.rds"))   

    ## work around for duplicate ISO A2 L2 names
    cases_per_district_samps <- merge_duplicate_district_data(cases_per_district_samps)
    pop_per_district <- merge_duplicate_district_data(pop_per_district)
    
    ## get the median annual incidence per district
    inc_per_district <- apply(cases_per_district_samps/pop_per_district, 2, median)
    cases_per_district <- apply(cases_per_district_samps, 2, median)

    ## rank districts by decreasing incidence
    district_rank_df <- tbl_df(data.frame(
        cntry_code = cntry,
        year = year,
        district = names(cases_per_district)[order(inc_per_district, decreasing=T)],
        incidence = inc_per_district[order(inc_per_district, decreasing=T)],
        cases = cases_per_district[order(inc_per_district, decreasing=T)],
        pop = pop_per_district[1,order(inc_per_district, decreasing=T)],
        index = order(inc_per_district, decreasing=T))) %>%
      dplyr::mutate(skip = ifelse(district %in% skip_ls | is.na(pop), TRUE, FALSE))
  return(district_rank_df)

}


#' @name identify_optimumDistrict_targets
#' @title identify_optimumDistrict_targets
#' @description Allocate vaccines according to the case optimized vaccine deployment strategy
#' @param country_ls
#' @param global_vacc_avail Total number of vaccines available to allocate, globally 
#' @param target_coverage_prop Vaccination coverage 
#' @param doses_per_person
#' @param rotate_every_nYears
#' @param year
#' @param out_wd
#' @param gen_out_wd
#' @return list with two dataframes, one for district-level (fullDistrict) targets and one for cell-level (partialDistrict) targets
identify_optimumDistrict_targets <- function(cntryCode_ls,
                              scenario,
                              global_vacc_avail,
                              target_coverage_prop,
                              doses_per_person = 2,
                              rotate_every_nYears = 3,
                              year,
                              partial_cover = FALSE,
                              partial_cover_method = "raster",
                              out_wd = "country_outputs/",
                              gen_out_wd = "country_targets/"
                              ){

    cat(sprintf("allocate vaccines with optimum district targeting by cases for %s: %s \n", year, paste(cntryCode_ls, collapse=", ")))
    
    inputs <- data.frame(country = cntryCode_ls, year = year) %>% 
        dplyr::mutate(country_year = paste(country, year, sep="--")) %>%
        dplyr::select(country_year) %>% unname %>% unlist
    all_district_ranks <- dplyr::bind_rows(mapply(FUN = rank_districts_by_country, country_year=inputs, scenario=scenario, rotate_every_nYears=rotate_every_nYears, in_wd = out_wd, out_wd = gen_out_wd, SIMPLIFY = FALSE)) %>%
        dplyr::filter(!skip)
    
    ## allocate vaccine to districts with highest cases
    all_district_ranks <- tbl_df(all_district_ranks[order(all_district_ranks$cases, decreasing = TRUE),]) %>%
      dplyr::mutate(id = seq_along(pop),
                    vacc_pop_coverage = pop * as.numeric(target_coverage_prop) * doses_per_person,
                    vacc_pop_cumsum = cumsum(vacc_pop_coverage),
                    vacc_used = ifelse(vacc_pop_cumsum <= global_vacc_avail, vacc_pop_coverage, 0),
                    fvp = vacc_used/doses_per_person,
                    fullDistrict_target = ifelse(vacc_used > 0, TRUE, FALSE))
    
    fullDistrict_target_names <- all_district_ranks$district[which(all_district_ranks$fullDistrict_target)]
    cat(paste0("Fully-targeted districts: ", paste(fullDistrict_target_names, collapse = ", "),"\n"))


    ## calculate number of vaccines remaining for district receiving partial allocation
    partialDistrict_vacc_pop_size <- global_vacc_avail - sum(all_district_ranks$vacc_used, na.rm = TRUE)
    ## identify highest ranked district that did not receive vaccine as the district that will receive a partial allocation
    ## NOTE: adding the criteria that pop needs to be greater than 0 to be targeted. This can occur with small polygons
    partialDistrict_target_id <- all_district_ranks$id[min(which(all_district_ranks$vacc_used==0 & all_district_ranks$pop>0))]
    partialDistrict_target_name <- all_district_ranks$district[partialDistrict_target_id]
    print(partialDistrict_vacc_pop_size)
    print(partialDistrict_target_id)
    cat(sprintf("Partially-targeted district: %s (%.1f vaccinated) \n ",
                partialDistrict_target_name, partialDistrict_vacc_pop_size))

    ## import cell-level rasters for the partial district country
    sublevel <- 2
    cntry_code <- all_district_ranks$cntry_code[partialDistrict_target_id]
    partialDistrict_target_poly_id <- all_district_ranks$index[partialDistrict_target_id] # index is the cntry-specific polygon ID
    cntry_code <- fix_country_name(cntry_code)
    lambda_per_cell <- R6raster$new(paste0(out_wd,cntry_code,"_lambda_per_cell.tif"))
    pop_per_cell <- R6raster$new(paste0(out_wd, cntry_code,"_",year, "_pop_per_cell.tif"))
    cases_per_cell <- R6raster$new(lambda_per_cell$raster * pop_per_cell$raster)
    cntry_shp <- get_country_sublevels(cntry_code, sublevel)
    
    wpRegion <- lookup_WorldPop_region(cntry_code)
    
    ## get rasters for partially targeted district with appropriate partial cover
    partialDistrict_lambda_rast <- extract_country_from_raster(lambda_per_cell, cntry_shp$geometry[partialDistrict_target_poly_id], partial_cover = FALSE, method = "raster")
    partialDistrict_pop_rast <- extract_country_from_raster(pop_per_cell, cntry_shp$geometry[partialDistrict_target_poly_id], partial_cover = partial_cover, method = partial_cover_method, trim_small = 0.05)
    partialDistrict_cases_rast <- R6raster$new(partialDistrict_lambda_rast$raster * partialDistrict_pop_rast$raster)

    ## allocate vaccine to the grid cells with greatest number of cases in partialDistrict at a given coverage level until there is no more
    partialDistrict_cases_summ_rast <- aggregate_raster_xlayers(partialDistrict_cases_rast, median, na.rm = TRUE)

    partialDistrict_lambda_summ_rast <- aggregate_raster_xlayers(partialDistrict_lambda_rast, median, na.rm = TRUE)

    ## some cells are NA, these ranks need to be last
    cell_rank <- data.frame(
        cntry_code = cntry_code,
        district = all_district_ranks$district[partialDistrict_target_id],
        year = year,
        index = order(values(partialDistrict_cases_summ_rast$raster), decreasing=T, na.last=TRUE),
        incidence = values(partialDistrict_lambda_summ_rast$raster)[order(values(partialDistrict_cases_summ_rast$raster), decreasing=T, na.last=TRUE)],
        cases = values(partialDistrict_cases_summ_rast$raster)[order(values(partialDistrict_cases_summ_rast$raster), decreasing=T, na.last=TRUE)],
        pop = values(partialDistrict_pop_rast$raster)[order(values(partialDistrict_cases_summ_rast$raster), decreasing=T, na.last=TRUE)])

    ## get all but the last cell to target here
    cell_rank <- cell_rank %>%
      dplyr::filter(!is.na(pop)) %>% ## added 6/3/18
      dplyr::mutate(vacc_pop_coverage = pop * target_coverage_prop * doses_per_person,
                  vacc_pop_cumsum = cumsum(vacc_pop_coverage),
                  vacc_used = ifelse(vacc_pop_cumsum <= partialDistrict_vacc_pop_size, vacc_pop_coverage, 0), 
                  fullCell_target = ifelse(vacc_used > 0, TRUE, FALSE))

    partialDistrict_target_cell_ids <- cell_rank[which(cell_rank$fullCell_target),]$index

    ## allocate last vaccines to highest-ranked untargeted cell, partially
    partialCell_vacc_pop_size <- partialDistrict_vacc_pop_size - sum(cell_rank$vacc_used, na.rm = TRUE)

    cat(paste0(" Full district vaccination: ", sum(all_district_ranks$vacc_used, na.rm = TRUE),"\n"))
    cat(paste0(" Partial district vaccination: ", partialDistrict_vacc_pop_size,"\n"))
    cat(paste0(" Total vaccination: ", sum(all_district_ranks$vacc_used, na.rm = TRUE) + partialDistrict_vacc_pop_size,"\n"))

    partialCell_target_cell_id <-cell_rank$index[min(which(cell_rank$vacc_used == 0), na.rm=T)] 

    ## workaround control flow for discrepancies between population at district and cell levels (should be fixed with 1km rasters - 6/4/18)
    if(!is.na(partialCell_target_cell_id)){
      cell_rank[which(cell_rank$index==partialCell_target_cell_id),]$vacc_used <- partialCell_vacc_pop_size
      partialCell_target_prop <- (partialCell_vacc_pop_size/doses_per_person)/cell_rank[which(cell_rank$index==partialCell_target_cell_id),]$pop

      print(paste("partial cell target ID", partialCell_target_cell_id))
      print(paste("vacc pop size", partialCell_vacc_pop_size))
      print(paste("target prop", partialCell_target_prop))
    } else{
      print(paste(partialCell_vacc_pop_size, "partial cell vaccines not allocated."))
    }

    cell_rank <- cell_rank %>%
      dplyr::mutate(fvp = vacc_used/doses_per_person,
                    prop_vaccinated = ifelse(pop==0, 0, fvp/pop))

    rm(lambda_per_cell, cases_per_cell, pop_per_cell)
    gc()

    if(!dir.exists(gen_out_wd)){
        dir.create(gen_out_wd)
    }

    ## Need these outputs to create vac_template_rast by country in a separate function
    write_csv(all_district_ranks, paste0(gen_out_wd, scenario, "_", year, "_districtTargets.csv"))
    write_csv(cell_rank, paste0(gen_out_wd, scenario, "_", year, "_cellTargets.csv"))

    return(list(district_targets = all_district_ranks, cell_targets = cell_rank))
}


#' @name identify_allDistrict_targets
#' @title identify_allDistrict_targets
#' @description Allocate vaccines according to the untargeted vaccine deployment strategy
#' @param country_ls
#' @param scenario
#' @param global_vacc_avail Total number of vaccines available to allocate, globally 
#' @param target_coverage_prop Vaccination coverage 
#' @param doses_per_person
#' @param rotate_every_nYears
#' @param year
#' @param out_wd
#' @param gen_out_wd
#' @return list with one dataframe, for district-level targets
identify_allDistrict_targets <- function(cntryCode_ls,
                              scenario,
                              global_vacc_avail,
                              target_coverage_prop, 
                              doses_per_person = 2,
                              rotate_every_nYears = 3, 
                              year,
                              out_wd = "country_outputs/",
                              gen_out_wd = "country_targets/"
                                  ){

    ## N.B. target_coverage_prop doesn't apply to this function
    cat(sprintf("allocate vaccines across all districts for %s: %s \n", year, paste(cntryCode_ls, collapse=", ")))
    
    inputs <- data.frame(country = cntryCode_ls, year = year) %>% 
        dplyr::mutate(country_year = paste(country, year, sep="--")) %>%
        dplyr::select(country_year) %>% unname %>% unlist
    all_district_noranks <- dplyr::bind_rows(mapply(FUN = rank_districts_by_country, country_year=inputs, scenario=scenario, rotate_every_nYears=rotate_every_nYears, in_wd = out_wd, out_wd = gen_out_wd, SIMPLIFY = FALSE)) ## nothing to skip for all targets

    ## calculate proportion of population that should be vaccinated based on vaccine availability
    global_vacc_prop <- global_vacc_avail/(sum(all_district_noranks$pop, na.rm=TRUE) * doses_per_person)

    all_district_noranks <- all_district_noranks %>%
        dplyr::mutate(id = seq_along(1:nrow(.)),
                      prop_vaccinated = global_vacc_prop,
                      vacc_used = pop * doses_per_person * global_vacc_prop,
                      fvp = vacc_used/doses_per_person)

    cat(paste0(" Total vaccination: ",sum(all_district_noranks$vacc_used, na.rm = TRUE),"\n"))
    cat(paste0(" Global pop vaccination proportion: ",global_vacc_prop,"\n"))

    if(!dir.exists(gen_out_wd)){
        dir.create(gen_out_wd)
    }

    write_csv(all_district_noranks, paste0(gen_out_wd, scenario, "_", year, "_districtTargets.csv"))

  return(list(district_targets = all_district_noranks))

}


#' @name rank_districts_by_country
#' @title rank_districts_by_country
#' @description Rank districts by total number of cases within a single country
#' @param country_year
#' @param scenario
#' @param rotate_every_nYears
#' @param in_wd
#' @param out_wd
#' @return dataframe with ordered districts and cases by country
rank_districts_by_country <- function(country_year, 
                                  scenario,
                                  rotate_every_nYears = 3,
                                  in_wd = "country_outputs/",
                                  out_wd = "country_targets/"
                                  ){

    country <- unlist(strsplit(country_year, split="--"))[1]
    year <- as.numeric(unlist(strsplit(country_year, split="--"))[2])

    ## get list of fully targeted districts in previous x years
    skipYears <- year-(1:(rotate_every_nYears-1))
    skip_ls <- skip_district_targets(skipYears, scenario, out_wd)

    cntry <- fix_country_name(country)
    ## import district-level data
    cases_per_district_samps <- readRDS(file = paste0(in_wd,cntry,"_",year,"_cases_per_district.rds"))
    pop_per_district <- readRDS(file = paste0(in_wd,cntry,"_",year,"_pop_per_district.rds"))   

    ## work around for duplicate ISO A2 L2 names
    cases_per_district_samps <- merge_duplicate_district_data(cases_per_district_samps)
    pop_per_district <- merge_duplicate_district_data(pop_per_district)

    ## get the median annual incidence per district
    inc_per_district <- apply(cases_per_district_samps/pop_per_district, 2, median)
    cases_per_district <- apply(cases_per_district_samps, 2, median)

    ## rank districts by decreasing cases
    district_rank_df <- data.frame(
        cntry_code = cntry,
        year = year,
        district = names(cases_per_district)[order(cases_per_district, decreasing=T)],
        incidence = inc_per_district[order(cases_per_district, decreasing=T)],
        cases = cases_per_district[order(cases_per_district, decreasing=T)],
        pop = pop_per_district[1,order(cases_per_district, decreasing=T)],
        index = order(cases_per_district, decreasing=T)) %>%
      dplyr::mutate(skip = ifelse(district %in% skip_ls | is.na(pop), TRUE, FALSE)) ## skip districts with missing pops
    return(district_rank_df)
}


#' @name identify_hrDistrict2_targets
#' @title identify_hrDistrict2_targets
#' @description Allocate vaccines to countries according to the case-logistics vaccination deployment strategy
#' @param cntryCode_ls
#' @param global_vacc_avail Total number of vaccines available to allocate, globally 
#' @param target_coverage_prop Vaccination coverage
#' @param doses_per_person
#' @param rotate_every_nYears
#' @param year
#' @param out_wd
#' @param gen_out_wd
#' @return list with two dataframes, one for district-level (fullDistrict) targets and one for cell-level (partialDistrict) targets
identify_hrDistrict2_targets <- function(cntryCode_ls,
                              scenario,
                              global_vacc_avail,
                              target_coverage_prop,
                              doses_per_person = 2,
                              rotate_every_nYears = 3,
                              year,
                              partial_cover = FALSE,
                              partial_cover_method = "raster",
                              out_wd = "country_outputs/",
                              gen_out_wd = "country_targets/"
                              ){

    cat(sprintf("allocate vaccines with 2-level high-risk district targeting for %s: %s \n", year, paste(cntryCode_ls, collapse=", ")))
    
    inputs <- data.frame(country = cntryCode_ls, year = year) %>% 
        dplyr::mutate(country_year = paste(country, year, sep="--")) %>%
        dplyr::select(country_year) %>% unname %>% unlist
    all_hrDistrict_ranks_orig <- pmap_dfr(.l = list(country_year=inputs, out_wd=out_wd), .f = get_hrDistricts2_by_country) 
    
    print(head(all_hrDistrict_ranks_orig)) 
    skipYears <- year-(1:(rotate_every_nYears-1))
    skip_ls <- skip_hrDistrict_targets2(skipYears, scenario, out_wd = gen_out_wd)

    ## ranked by total high-risk population across years
    countries_ranked_by_hrDistrict1 <- all_hrDistrict_ranks_orig %>%
      dplyr::filter(hr_district & hr_level == 1) %>%
      dplyr::mutate(hr_incidence = hr_pop/pop) %>%
      group_by(cntry_code) %>%
      summarise(hr_level = first(hr_level), hr_pop = sum(hr_pop, na.rm = TRUE), pop = first(pop), cases = sum(cases, na.rm = TRUE)) %>%
      arrange(desc(cases)) %>%
      dplyr::mutate(cntry_rank = seq_along(cases)) %>%
      ungroup %>% 
      dplyr::select(cntry_code, cntry_rank)
    hr1_countries <- unique(countries_ranked_by_hrDistrict1$cntry_code)
    
    ## countries with only 2nd level high-risk populations
    countries_ranked_by_hrDistrict2 <- all_hrDistrict_ranks_orig %>%
      dplyr::filter(hr_district & hr_level == 2 & !(cntry_code %in% hr1_countries)) %>%
      dplyr::mutate(hr_incidence = hr_pop/pop) %>%
      group_by(cntry_code) %>%
      summarise(hr_level = first(hr_level), hr_pop = sum(hr_pop, na.rm = TRUE), pop = first(pop), cases = sum(cases, na.rm = TRUE)) %>%
      arrange(desc(cases)) %>%
      dplyr::mutate(cntry_rank = seq_along(cases)+nrow(countries_ranked_by_hrDistrict1)) %>%
      ungroup %>% 
      dplyr::select(cntry_code, cntry_rank)
    hr2_countries <- unique(countries_ranked_by_hrDistrict2$cntry_code)
  
    ## countries with only 3rd level high-risk populations
    countries_ranked_by_hrDistrict3 <- all_hrDistrict_ranks_orig %>%
      dplyr::filter(hr_district & hr_level == 3 & !(cntry_code %in% hr1_countries) & !(cntry_code %in% hr2_countries)) %>%
      dplyr::mutate(hr_incidence = hr_pop/pop) %>%
      group_by(cntry_code) %>%
      summarise(hr_level = first(hr_level), hr_pop = sum(hr_pop, na.rm = TRUE), pop = first(pop), cases = sum(cases, na.rm = TRUE)) %>%
      arrange(desc(cases)) %>%
      dplyr::mutate(cntry_rank = seq_along(cases)+nrow(countries_ranked_by_hrDistrict1)+nrow(countries_ranked_by_hrDistrict2)) %>%
      ungroup %>% 
      dplyr::select(cntry_code, cntry_rank)
    hr3_countries <- unique(countries_ranked_by_hrDistrict3$cntry_code)
    
    ## countries with only 4th level high-risk populations
    countries_ranked_by_hrDistrict4 <- all_hrDistrict_ranks_orig %>%
      dplyr::filter(hr_district & hr_level == 4 & !(cntry_code %in% hr1_countries) & !(cntry_code %in% hr2_countries) & !(cntry_code %in% hr3_countries)) %>%
      dplyr::mutate(hr_incidence = hr_pop/pop) %>%
      group_by(cntry_code) %>%
      summarise(hr_level = first(hr_level), hr_pop = sum(hr_pop, na.rm = TRUE), pop = first(pop), cases = sum(cases, na.rm = TRUE)) %>%
      arrange(desc(cases)) %>%
      dplyr::mutate(cntry_rank = seq_along(cases)+nrow(countries_ranked_by_hrDistrict1)+nrow(countries_ranked_by_hrDistrict2)+nrow(countries_ranked_by_hrDistrict3)) %>%
      ungroup %>% 
      dplyr::select(cntry_code, cntry_rank)
   
    ## bind 2 levels of high-risk pops
    countries_ranked_by_hrDistrict <- bind_rows(countries_ranked_by_hrDistrict1, countries_ranked_by_hrDistrict2, countries_ranked_by_hrDistrict3, countries_ranked_by_hrDistrict4)
    avail_hrDistrict_ranks <- full_join(all_hrDistrict_ranks_orig, countries_ranked_by_hrDistrict, by = c("cntry_code")) %>%
      dplyr::filter(hr_district) %>%
      dplyr::mutate(skip = ifelse(district %in% skip_ls | is.na(pop), TRUE, FALSE)) %>%
      dplyr::filter(!skip)
    print(countries_ranked_by_hrDistrict)
    print(avail_hrDistrict_ranks)
    print("countries ranked by hrDistrict2 tiers")

    ##############################################################
    #### skip targeting if there are no available hrDistricts ####
    if (nrow(avail_hrDistrict_ranks) > 0){

      all_hrDistrict_ranks <- avail_hrDistrict_ranks %>%
        arrange(hr_level, cntry_rank, desc(cases)) %>%
        dplyr::mutate(id = seq_along(cases),
                      vacc_pop_coverage = pop * as.numeric(target_coverage_prop) * doses_per_person,
                      vacc_pop_cumsum = cumsum(vacc_pop_coverage),
                      vacc_used = ifelse(vacc_pop_cumsum <= global_vacc_avail, vacc_pop_coverage, 0),
                      fvp = vacc_used/doses_per_person,
                      fullDistrict_target = ifelse(vacc_used > 0, TRUE, FALSE))

      fullDistrict_target_names <- all_hrDistrict_ranks$district[which(all_hrDistrict_ranks$fullDistrict_target)]
      cat(paste0("Fully-targeted districts: ", paste(fullDistrict_target_names, collapse = ", "),"\n"))

      ## calculate number of vaccines remaining for district receiving partial allocation
      partialDistrict_vacc_pop_size <- global_vacc_avail - sum(all_hrDistrict_ranks$vacc_used, na.rm = TRUE)
      ## identify highest ranked district that did not receive vaccine as the district that will receive a partial allocation
      ## NOTE: adding the criteria that pop needs to be greater than 0 to be targeted. This can occur with small polygons

      ## control flow for partial district targets
      if(sum(all_hrDistrict_ranks$fullDistrict_target, na.rm=TRUE) < nrow(all_hrDistrict_ranks)){

        partialDistrict_target_id <- all_hrDistrict_ranks$id[min(which(all_hrDistrict_ranks$vacc_used==0 & all_hrDistrict_ranks$pop>0))]
        partialDistrict_target_name <- all_hrDistrict_ranks$district[partialDistrict_target_id]
        print(partialDistrict_vacc_pop_size)
        print(partialDistrict_target_id)
        cat(sprintf("Partially-targeted district: %s (%.1f vaccinated) \n ",
                    partialDistrict_target_name, partialDistrict_vacc_pop_size))

        ## import cell-level rasters for the partial district country
        sublevel <- 2
        cntry_code <- all_hrDistrict_ranks$cntry_code[partialDistrict_target_id]
        partialDistrict_target_poly_id <- all_hrDistrict_ranks$index[partialDistrict_target_id] # index is the cntry-specific polygon ID
        cntry_code <- fix_country_name(cntry_code)
        lambda_per_cell <- R6raster$new(paste0(out_wd,cntry_code,"_lambda_per_cell.tif"))
        pop_per_cell <- R6raster$new(paste0(out_wd, cntry_code,"_",year, "_pop_per_cell.tif"))
        cases_per_cell <- R6raster$new(lambda_per_cell$raster * pop_per_cell$raster)
        cntry_shp <- get_country_sublevels(cntry_code, sublevel)

        ## get rasters for partially targeted district
        partialDistrict_lambda_rast <- extract_country_from_raster(lambda_per_cell, cntry_shp$geometry[partialDistrict_target_poly_id], partial_cover = FALSE, method = "raster")
        partialDistrict_pop_rast <- extract_country_from_raster(pop_per_cell, cntry_shp$geometry[partialDistrict_target_poly_id], partial_cover = partial_cover, method = partial_cover_method, trim_small = 0.05)
        partialDistrict_cases_rast <- R6raster$new(partialDistrict_lambda_rast$raster * partialDistrict_pop_rast$raster)

        ## allocate vaccine to the grid cells with greatest number of cases in partialDistrict at a given coverage level until there is no more
        partialDistrict_cases_summ_rast <- aggregate_raster_xlayers(partialDistrict_cases_rast, median, na.rm = TRUE)
       
        partialDistrict_lambda_summ_rast <- aggregate_raster_xlayers(partialDistrict_lambda_rast, median, na.rm = TRUE)

        ## some cells are NA, these ranks need to be last
        cell_rank <- data.frame(
            cntry_code = cntry_code,
            district = all_hrDistrict_ranks$district[partialDistrict_target_id],
            year = year,
            index = order(values(partialDistrict_cases_summ_rast$raster), decreasing=T, na.last=TRUE),
            incidence = values(partialDistrict_lambda_summ_rast$raster)[order(values(partialDistrict_cases_summ_rast$raster), decreasing=T, na.last=TRUE)],
            cases = values(partialDistrict_cases_summ_rast$raster)[order(values(partialDistrict_cases_summ_rast$raster), decreasing=T, na.last=TRUE)],
            pop = values(partialDistrict_pop_rast$raster)[order(values(partialDistrict_cases_summ_rast$raster), decreasing=T, na.last=TRUE)])

        ## get all but the last cell to target here
        cell_rank <- cell_rank %>%
        dplyr::mutate(vacc_pop_coverage = pop * as.numeric(target_coverage_prop) * doses_per_person,
                      vacc_pop_cumsum = cumsum(vacc_pop_coverage),
                      vacc_used = ifelse(vacc_pop_cumsum <= partialDistrict_vacc_pop_size, vacc_pop_coverage, 0),
                      fullCell_target = ifelse(vacc_used > 0, TRUE, FALSE)) %>%
        dplyr::mutate(fvp = vacc_used/doses_per_person,
                      prop_vaccinated = ifelse(pop==0, 0, fvp/pop))

        partialDistrict_target_cell_ids <- cell_rank[which(cell_rank$fullCell_target),]$index

        ## allocate last vaccines to highest-ranked untargeted cell, partially
        partialCell_vacc_pop_size <- partialDistrict_vacc_pop_size - sum(cell_rank$vacc_used, na.rm = TRUE)

        cat(paste0(" Full district vaccination: ", sum(all_hrDistrict_ranks$vacc_used, na.rm = TRUE),"\n"))
        cat(paste0(" Partial district vaccination: ", partialDistrict_vacc_pop_size,"\n"))
        cat(paste0(" Total vaccination: ", sum(all_hrDistrict_ranks$vacc_used, na.rm = TRUE) +
                                          partialDistrict_vacc_pop_size,"\n"))

        ## control flow for partial cell targeting
        if ((sum(cell_rank$fullCell_target, na.rm=TRUE)) < nrow(cell_rank %>% dplyr::filter(!is.na(pop)))){

          partialCell_target_cell_id <-cell_rank$index[min(which(cell_rank$vacc_used == 0 & cell_rank$pop>0), na.rm=T)]

          if(!is.na(partialCell_target_cell_id)){
            cell_rank[which(cell_rank$index==partialCell_target_cell_id),]$vacc_used <- partialCell_vacc_pop_size
            partialCell_target_prop <- (partialCell_vacc_pop_size/doses_per_person)/cell_rank[which(cell_rank$index==partialCell_target_cell_id),]$pop

            print(paste("partial cell target ID", partialCell_target_cell_id))
            print(paste("vacc pop size", partialCell_vacc_pop_size))
            print(paste("target prop", partialCell_target_prop))
          } else{
            print(paste(partialCell_vacc_pop_size, "partial cell vaccines not allocated."))
          }
          
          ## recalculate fvp and prop_vaccinated if there is partial cell targeting
          cell_rank <- cell_rank %>%
            dplyr::mutate(fvp = vacc_used/doses_per_person,
                          prop_vaccinated = ifelse(pop==0, 0, fvp/pop))
        
        } ## pc targeting

        rm(lambda_per_cell, cases_per_cell, pop_per_cell)
        gc()
      
      } ## pd targeting

      if(!dir.exists(gen_out_wd)){
        dir.create(gen_out_wd)
      }

      ## Need these outputs to create vac_template_rast by country in a separate function
      skipped_districts <- full_join(all_hrDistrict_ranks_orig, countries_ranked_by_hrDistrict, by = c("cntry_code")) %>%
        dplyr::filter(hr_district) %>%
        dplyr::mutate(skip = ifelse(district %in% skip_ls | is.na(pop), TRUE, FALSE)) %>%
        dplyr::filter(skip) %>%
        dplyr::mutate(id = NA,
                      vacc_pop_coverage = NA,
                      vacc_pop_cumsum = NA,
                      vacc_used = NA,
                      fvp = NA,
                      fullDistrict_target = NA)
      all_hrDistrict_ranks2 <- bind_rows(all_hrDistrict_ranks, skipped_districts)
      write_csv(all_hrDistrict_ranks2, paste0(gen_out_wd, scenario, "_", year, "_districtTargets.csv"))
      
      if(exists("cell_rank")){
        write_csv(cell_rank, paste0(gen_out_wd, scenario, "_", year, "_cellTargets.csv"))
        } else{
        cell_rank <- data.frame()
        }
      
      return(list(district_targets = all_hrDistrict_ranks, cell_targets = cell_rank))

    }  else { ## no hrTargets were found
      print("No high-risk targets are available")
      return(list(district_targets = data.frame(), cell_targets = data.frame()))
    }
    
}


#' @name get_hrDistricts2_by_country
#' @title get_hrDistricts2_by_country
#' @description Rank all high risk districts in a country by cases, where high risk districts are ones where the incidence rate is greater than 1/1000 for at least 100,000 residents or 10% of the district's population
#' @param country_year
#' @param hrThreshold1
#' @param hrThreshold2
#' @param hrThreshold3
#' @param hrThreshold4
#' @param partial_cover
#' @param partial_cover_method
#' @param out_wd
#' @return dataframe with ordered cells and incidence by country
get_hrDistricts2_by_country <- function(country_year,
                                  hrThreshold1 = 1/1000,
                                  hrThreshold2 = 1/5000,
                                  hrThreshold3 = 1/10000,
                                  hrThreshold4 = 1/100000,
                                  partial_cover = FALSE,
                                  partial_cover_method = "raster",
                                  out_wd= "manuscripts/GAVI Impact Estimation/country_outputs_ms/"
                                  ){

    country <- unlist(strsplit(country_year, split="--"))[1]
    year <- unlist(strsplit(country_year, split="--"))[2]
    cat(sprintf("ranking %s cells in %s \n", country, year))

    cntry <- fix_country_name(country)

    ## import pre-saved data
    lambda_per_cell <- R6raster$new(paste0(out_wd,cntry,"_lambda_per_cell.tif"))
    pop_per_district <- readRDS(paste0(out_wd,cntry,"_",year,"_pop_per_district.rds"))
    cases_per_district_samps <- readRDS(paste0(out_wd,cntry,"_",year,"_cases_per_district.rds"))
  
    print(colnames(cases_per_district_samps)[duplicated(colnames(cases_per_district_samps))])

    ## work around for duplicate ISO A2 L2 names
    cases_per_district_samps <- merge_duplicate_district_data(cases_per_district_samps)
    pop_per_district <- merge_duplicate_district_data(pop_per_district)

    inc_per_district <- apply(cases_per_district_samps/pop_per_district, 2, median)
    print("presaved data imported")
    who_region <- lookup_WorldPop_region(cntry)
    
    ## identify high risk districts
    ## summarize lambdas across samples
    cases_per_district_summ <- apply(cases_per_district_samps, 2, median, na.rm=TRUE)
    print("cases per district median")
    
    ## import existing high-risk population files if available
    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop_per_cell.tif"))){
      lambda_summ_rast <- aggregate_raster_xlayers(lambda_per_cell, median, na.rm = TRUE)
    print("lambda per cell, aggregated")
    ## import new pop version where partial_cover = FALSE
      pop_per_cell_npc <- load_population_estimates(region = who_region, lambda_per_cell, as.numeric(year))
    print("pop per cell npc, loaded")
      ## filter cells with high incidence (>1/1000)
      hr_lambda_summ_rast <- R6raster$new(raster::calc(lambda_summ_rast$raster, fun=function(x){ifelse(x > hrThreshold1, x, NA)}))
    print("hr lambda summ threshold calculated")
      ## get population for high risk cells (grab original data from non-partial cover version)
      hr_pop_per_cell <- R6raster$new(raster::mask(pop_per_cell_npc$raster, hr_lambda_summ_rast$raster, maskvalue=NA))
    print("hr pop per cell masked")

      writeRaster(hr_pop_per_cell$raster, paste0(out_wd,cntry,"_",year,"_hrpop_per_cell.tif"))
      
      rm(pop_per_cell_npc, hr_lambda_summ_rast)
      gc()

    } else {
      hr_pop_per_cell <- R6raster$new(paste0(out_wd,cntry,"_",year,"_hrpop_per_cell.tif"))
    print("hr pop per cell imported")
    }
   
    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop_per_district.rds"))){
      ## summarize hr pop by district
      hr_pop_per_district <- apply_to_all_sublevels(cntry, 
                                            hr_pop_per_cell, 
                                            ISO_level = 2, 
                                            partial_cover = partial_cover,
                                            method = partial_cover_method,
                                           trim_small = 0.05,
                                            fun=function(x){if(!all(is.na(x))){
                                              rc <- sum(x, na.rm=TRUE)
                                            } else {
                                              rc <- NA
                                            }})
      saveRDS(hr_pop_per_district, paste0(out_wd,cntry,"_",year,"_hrpop_per_district.rds"))
  print("hr pop per district saved")
    } else{
      hr_pop_per_district <- readRDS(paste0(out_wd,cntry,"_",year,"_hrpop_per_district.rds"))
  print("hr pop per district imported")
    }

  ## work around for duplicate ISO A2 L2 names
  hr_pop_per_district <- merge_duplicate_district_data(hr_pop_per_district)

    ## import second-level high-risk population files if available
    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop5k_per_cell.tif"))){
      lambda_summ_rast <- aggregate_raster_xlayers(lambda_per_cell, median, na.rm = TRUE)
    print("lambda per cell, aggregated")
    ## import new pop version where partial_cover = FALSE
      pop_per_cell_npc <- load_population_estimates(region = who_region, lambda_per_cell, as.numeric(year))
    print("pop per cell npc, loaded")
      ## filter cells with high incidence (>1/1000)
      hr2_lambda_summ_rast <- R6raster$new(raster::calc(lambda_summ_rast$raster, fun=function(x){ifelse(x > hrThreshold2, x, NA)}))
    print("hr2 lambda summ threshold calculated")
      ## get population for high risk cells (grab original data from non-partial cover version)
      hr2_pop_per_cell <- R6raster$new(raster::mask(pop_per_cell_npc$raster, hr2_lambda_summ_rast$raster, maskvalue=NA))
    print("hr2 pop per cell masked")

      writeRaster(hr2_pop_per_cell$raster, paste0(out_wd,cntry,"_",year,"_hrpop5k_per_cell.tif"))
      
      rm(pop_per_cell_npc, hr2_lambda_summ_rast)
      gc()
    } else{
      hr2_pop_per_cell <- R6raster$new(paste0(out_wd,cntry,"_",year,"_hrpop5k_per_cell.tif"))
  print("hr2 pop per cell imported")
    }

    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop5k_per_district.rds"))){
      ## summarize hr pop by district
      hr2_pop_per_district <- apply_to_all_sublevels(cntry, 
                                            hr2_pop_per_cell, 
                                            ISO_level = 2, 
                                            partial_cover = partial_cover,
                                            method = partial_cover_method,
                                           trim_small = 0.05,
                                            fun=function(x){if(!all(is.na(x))){
                                              rc <- sum(x, na.rm=TRUE)
                                            } else {
                                              rc <- NA
                                            }})
      saveRDS(hr2_pop_per_district, paste0(out_wd,cntry,"_",year,"_hrpop5k_per_district.rds"))
  print("hr2 pop per district saved")
    } else{
      hr2_pop_per_district <- readRDS(paste0(out_wd,cntry,"_",year,"_hrpop5k_per_district.rds"))
  print("hr2 pop per district imported")
    }

    ## work around for duplicate ISO A2 L2 names
    hr2_pop_per_district <- merge_duplicate_district_data(hr2_pop_per_district)

    ## import third-level high-risk population files if available
    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop10k_per_cell.tif"))){
      lambda_summ_rast <- aggregate_raster_xlayers(lambda_per_cell, median, na.rm = TRUE)
    print("lambda per cell, aggregated")
    ## import new pop version where partial_cover = FALSE
      pop_per_cell_npc <- load_population_estimates(region = who_region, lambda_per_cell, as.numeric(year))
    print("pop per cell npc, loaded")
      ## filter cells with high incidence (>1/10000)
      hr3_lambda_summ_rast <- R6raster$new(raster::calc(lambda_summ_rast$raster, fun=function(x){ifelse(x > hrThreshold3, x, NA)}))
    print("hr3 lambda summ threshold calculated")
      ## get population for high risk cells (grab original data from non-partial cover version)
      hr3_pop_per_cell <- R6raster$new(raster::mask(pop_per_cell_npc$raster, hr3_lambda_summ_rast$raster, maskvalue=NA))
    print("hr3 pop per cell masked")

      writeRaster(hr3_pop_per_cell$raster, paste0(out_wd,cntry,"_",year,"_hrpop10k_per_cell.tif"))
      
      rm(pop_per_cell_npc, hr3_lambda_summ_rast)
      gc()
    } else{
      hr3_pop_per_cell <- R6raster$new(paste0(out_wd,cntry,"_",year,"_hrpop10k_per_cell.tif"))
  print("hr3 pop per cell imported")
    }

    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop10k_per_district.rds"))){
      ## summarize hr pop by district
      hr3_pop_per_district <- apply_to_all_sublevels(cntry, 
                                            hr3_pop_per_cell, 
                                            ISO_level = 2, 
                                            partial_cover = partial_cover,
                                            method = partial_cover_method,
                                           trim_small = 0.05,
                                            fun=function(x){if(!all(is.na(x))){
                                              rc <- sum(x, na.rm=TRUE)
                                            } else {
                                              rc <- NA
                                            }})
      saveRDS(hr3_pop_per_district, paste0(out_wd,cntry,"_",year,"_hrpop10k_per_district.rds"))
  print("hr3 pop per district saved")
    } else{
      hr3_pop_per_district <- readRDS(paste0(out_wd,cntry,"_",year,"_hrpop10k_per_district.rds"))
  print("hr3 pop per district imported")
    }

    ## work around for duplicate ISO A2 L2 names
    hr3_pop_per_district <- merge_duplicate_district_data(hr3_pop_per_district)

    ## import fourth-level high-risk population files if available
    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop100k_per_cell.tif"))){
      lambda_summ_rast <- aggregate_raster_xlayers(lambda_per_cell, median, na.rm = TRUE)
    print("lambda per cell, aggregated")
    ## import new pop version where partial_cover = FALSE
      pop_per_cell_npc <- load_population_estimates(region = who_region, lambda_per_cell, as.numeric(year))
    print("pop per cell npc, loaded")
      ## filter cells with high incidence (>1/10000)
      hr4_lambda_summ_rast <- R6raster$new(raster::calc(lambda_summ_rast$raster, fun=function(x){ifelse(x > hrThreshold4, x, NA)}))
    print("hr4 lambda summ threshold calculated")
      ## get population for high risk cells (grab original data from non-partial cover version)
      hr4_pop_per_cell <- R6raster$new(raster::mask(pop_per_cell_npc$raster, hr4_lambda_summ_rast$raster, maskvalue=NA))
    print("hr4 pop per cell masked")

      writeRaster(hr4_pop_per_cell$raster, paste0(out_wd,cntry,"_",year,"_hrpop100k_per_cell.tif"))
      
      rm(pop_per_cell_npc, hr4_lambda_summ_rast)
      gc()
    } else{
      hr4_pop_per_cell <- R6raster$new(paste0(out_wd,cntry,"_",year,"_hrpop100k_per_cell.tif"))
  print("hr4 pop per cell imported")
    }

    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop100k_per_district.rds"))){
      ## summarize hr pop by district
      hr4_pop_per_district <- apply_to_all_sublevels(cntry, 
                                            hr4_pop_per_cell, 
                                            ISO_level = 2, 
                                            partial_cover = partial_cover,
                                            method = partial_cover_method,
                                           trim_small = 0.05,
                                            fun=function(x){if(!all(is.na(x))){
                                              rc <- sum(x, na.rm=TRUE)
                                            } else {
                                              rc <- NA
                                            }})
      saveRDS(hr4_pop_per_district, paste0(out_wd,cntry,"_",year,"_hrpop100k_per_district.rds"))
  print("hr4 pop per district saved")
    } else{
      hr4_pop_per_district <- readRDS(paste0(out_wd,cntry,"_",year,"_hrpop100k_per_district.rds"))
  print("hr4 pop per district imported")
    }

    ## work around for duplicate ISO A2 L2 names
    hr4_pop_per_district <- merge_duplicate_district_data(hr4_pop_per_district)
    
    ## highest risk districts
    hr1_df <- tbl_df(data.frame(cntry_code = cntry, 
                                  year = year, 
                                  district = colnames(pop_per_district)[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  pop = pop_per_district[1, order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  hr_pop = hr_pop_per_district[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  cases = cases_per_district_summ[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  incidence = inc_per_district[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  index = order(cases_per_district_summ, decreasing=T, na.last=TRUE))) %>%
      dplyr::mutate(hr_prop = ifelse(pop > 0, hr_pop/pop, 0),
                    hr_district = ifelse(hr_prop > 0.1 | hr_pop > 100000, TRUE, FALSE),
                    hr_level = 1)
    hr_dnames <- hr1_df %>% dplyr::filter(!is.na(hr_pop)) %>% dplyr::distinct(district) %>% unlist

    ## second-level high risk districts
    hr2_df <- tbl_df(data.frame(cntry_code = cntry, 
                                  year = year, 
                                  district = colnames(pop_per_district)[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  pop = pop_per_district[1, order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  hr_pop = hr2_pop_per_district[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  cases = cases_per_district_summ[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  incidence = inc_per_district[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  index = order(cases_per_district_summ, decreasing=T, na.last=TRUE))) %>%
      dplyr::mutate(hr_prop = ifelse(pop > 0, hr_pop/pop, 0),
                    hr_district = ifelse(hr_prop > 0.1 | hr_pop > 100000, TRUE, FALSE),
                    hr_level = 2) %>%
      dplyr::filter(!(district %in% hr_dnames)) ## rm hrDistrict level 1
    hr_dnames2 <- hr2_df %>% dplyr::filter(!is.na(hr_pop)) %>% dplyr::distinct(district) %>% unlist


    ## third-level high risk districts
    hr3_df <- tbl_df(data.frame(cntry_code = cntry, 
                                  year = year, 
                                  district = colnames(pop_per_district)[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  pop = pop_per_district[1, order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  hr_pop = hr3_pop_per_district[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  cases = cases_per_district_summ[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  incidence = inc_per_district[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  index = order(cases_per_district_summ, decreasing=T, na.last=TRUE))) %>%
      dplyr::mutate(hr_prop = ifelse(pop > 0, hr_pop/pop, 0),
                    hr_district = ifelse(hr_prop > 0.1 | hr_pop > 100000, TRUE, FALSE),
                    hr_level = 3) %>%
      dplyr::filter(!(district %in% hr_dnames) & !(district %in% hr_dnames2)) ## rm hrDistrict level 1 & 2
    hr_dnames3 <- hr3_df %>% dplyr::filter(!is.na(hr_pop)) %>% dplyr::distinct(district) %>% unlist

    ## fourth-level high risk districts
    hr4_df <- tbl_df(data.frame(cntry_code = cntry, 
                                  year = year, 
                                  district = colnames(pop_per_district)[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  pop = pop_per_district[1, order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  hr_pop = hr4_pop_per_district[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  cases = cases_per_district_summ[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  incidence = inc_per_district[order(cases_per_district_summ, decreasing=T, na.last=TRUE)], 
                                  index = order(cases_per_district_summ, decreasing=T, na.last=TRUE))) %>%
      dplyr::mutate(hr_prop = ifelse(pop > 0, hr_pop/pop, 0),
                    hr_district = ifelse(hr_prop > 0.1 | hr_pop > 100000, TRUE, FALSE),
                    hr_level = 4) %>%
      dplyr::filter(!(district %in% hr_dnames) & !(district %in% hr_dnames2) & !(district %in% hr_dnames3)) ## rm hrDistrict level 1 & 2 & 3

    hr_df <- bind_rows(hr1_df, hr2_df, hr3_df, hr4_df)
    
  rm(lambda_per_cell, hr_pop_per_cell)
  gc()

  return(hr_df)
}


#' @name identify_hrDistrict2Incid_targets
#' @title identify_hrDistrict2Incid_targets
#' @description Allocate vaccines to countries according to the rate-logistics vaccination deployment strategy
#' @param cntryCode_ls
#' @param global_vacc_avail Total number of vaccines available to allocate, globally 
#' @param target_coverage_prop Vaccination coverage 
#' @param doses_per_person
#' @param rotate_every_nYears
#' @param year
#' @param out_wd
#' @param gen_out_wd
#' @return list with two dataframes, one for district-level (fullDistrict) targets and one for cell-level (partialDistrict) targets
identify_hrDistrict2Incid_targets <- function(cntryCode_ls,
                              scenario,
                              global_vacc_avail,
                              target_coverage_prop,
                              doses_per_person = 2,
                              rotate_every_nYears = 3,
                              year,
                              partial_cover = FALSE,
                              partial_cover_method = "raster",
                              out_wd = "country_outputs/",
                              gen_out_wd = "country_targets/"
                              ){

    cat(sprintf("allocate vaccines with 2-level high-risk district targeting for %s: %s \n", year, paste(cntryCode_ls, collapse=", ")))
    
    inputs <- data.frame(country = cntryCode_ls, year = year) %>% 
        dplyr::mutate(country_year = paste(country, year, sep="--")) %>%
        dplyr::select(country_year) %>% unname %>% unlist
    all_hrDistrict_ranks_orig <- pmap_dfr(.l = list(country_year=inputs, out_wd=out_wd), .f = get_hrDistricts2Incid_by_country) 
    
      print(head(all_hrDistrict_ranks_orig)) 
    skipYears <- year-(1:(rotate_every_nYears-1))
    ## skip_hrDistrict_targets1 lists hrDistricts only for countries that were fully targeted in previous campaigns
    ## skip_hrDistrict_targets2 lists all hrDistricts that were fully targeted in previous campaigns, regardless of country status 
    skip_ls <- skip_hrDistrict_targets2(skipYears, scenario, out_wd = gen_out_wd)

    ## ranked by total high-risk population across years
    countries_ranked_by_hrDistrict1 <- all_hrDistrict_ranks_orig %>%
      dplyr::filter(hr_district & hr_level == 1) %>%
      dplyr::mutate(hr_incidence = hr_pop/pop) %>%
      group_by(cntry_code) %>%
      summarise(hr_level = first(hr_level), hr_pop = sum(hr_pop, na.rm = TRUE), pop = first(pop), cases = sum(cases, na.rm = TRUE)) %>%
      arrange(desc(cases)) %>%
      dplyr::mutate(cntry_rank = seq_along(cases)) %>%
      ungroup %>% 
      dplyr::select(cntry_code, cntry_rank)
    hr1_countries <- unique(countries_ranked_by_hrDistrict1$cntry_code)
      
    ## countries with only 2nd level high-risk populations
    countries_ranked_by_hrDistrict2 <- all_hrDistrict_ranks_orig %>%
      dplyr::filter(hr_district & hr_level == 2 & !(cntry_code %in% hr1_countries)) %>%
      dplyr::mutate(hr_incidence = hr_pop/pop) %>%
      group_by(cntry_code) %>%
      summarise(hr_level = first(hr_level), hr_pop = sum(hr_pop, na.rm = TRUE), pop = first(pop), cases = sum(cases, na.rm = TRUE)) %>%
      arrange(desc(cases)) %>%
      dplyr::mutate(cntry_rank = seq_along(cases)+nrow(countries_ranked_by_hrDistrict1)) %>%
      ungroup %>% 
      dplyr::select(cntry_code, cntry_rank)
    hr2_countries <- unique(countries_ranked_by_hrDistrict2$cntry_code)
    
    ## countries with only 3rd level high-risk populations
    countries_ranked_by_hrDistrict3 <- all_hrDistrict_ranks_orig %>%
      dplyr::filter(hr_district & hr_level == 3 & !(cntry_code %in% hr1_countries) & !(cntry_code %in% hr2_countries)) %>%
      dplyr::mutate(hr_incidence = hr_pop/pop) %>%
      group_by(cntry_code) %>%
      summarise(hr_level = first(hr_level), hr_pop = sum(hr_pop, na.rm = TRUE), pop = first(pop), cases = sum(cases, na.rm = TRUE)) %>%
      arrange(desc(cases)) %>%
      dplyr::mutate(cntry_rank = seq_along(cases)+nrow(countries_ranked_by_hrDistrict1)+nrow(countries_ranked_by_hrDistrict2)) %>%
      ungroup %>% 
      dplyr::select(cntry_code, cntry_rank)
    hr3_countries <- unique(countries_ranked_by_hrDistrict3$cntry_code)
      
    ## countries with only 4th level high-risk populations
    countries_ranked_by_hrDistrict4 <- all_hrDistrict_ranks_orig %>%
      dplyr::filter(hr_district & hr_level == 4 & !(cntry_code %in% hr1_countries) & !(cntry_code %in% hr2_countries) & !(cntry_code %in% hr3_countries)) %>%
      dplyr::mutate(hr_incidence = hr_pop/pop) %>%
      group_by(cntry_code) %>%
      summarise(hr_level = first(hr_level), hr_pop = sum(hr_pop, na.rm = TRUE), pop = first(pop), cases = sum(cases, na.rm = TRUE)) %>%
      arrange(desc(cases)) %>%
      dplyr::mutate(cntry_rank = seq_along(cases)+nrow(countries_ranked_by_hrDistrict1)+nrow(countries_ranked_by_hrDistrict2)+nrow(countries_ranked_by_hrDistrict3)) %>%
      ungroup %>% 
      dplyr::select(cntry_code, cntry_rank)
     
    ## bind 2 levels of high-risk pops
    countries_ranked_by_hrDistrict <- bind_rows(countries_ranked_by_hrDistrict1, countries_ranked_by_hrDistrict2, countries_ranked_by_hrDistrict3, countries_ranked_by_hrDistrict4)
    avail_hrDistrict_ranks <- full_join(all_hrDistrict_ranks_orig, countries_ranked_by_hrDistrict, by = c("cntry_code")) %>%
      dplyr::filter(hr_district) %>%
      dplyr::mutate(skip = ifelse(district %in% skip_ls | is.na(pop), TRUE, FALSE)) %>%
      dplyr::filter(!skip)
      print(countries_ranked_by_hrDistrict)
      print(avail_hrDistrict_ranks)
    print("countries ranked by hrDistrict2Incid tiers")

    ##############################################################
    #### skip targeting if there are no available hrDistricts ####
    if (nrow(avail_hrDistrict_ranks) > 0){

      all_hrDistrict_ranks <- avail_hrDistrict_ranks %>%
        arrange(hr_level, cntry_rank, desc(incidence)) %>%
        dplyr::mutate(id = seq_along(incidence),
                      vacc_pop_coverage = pop * as.numeric(target_coverage_prop) * doses_per_person,
                      vacc_pop_cumsum = cumsum(vacc_pop_coverage),
                      vacc_used = ifelse(vacc_pop_cumsum <= global_vacc_avail, vacc_pop_coverage, 0),
                      fvp = vacc_used/doses_per_person,
                      fullDistrict_target = ifelse(vacc_used > 0, TRUE, FALSE))

      fullDistrict_target_names <- all_hrDistrict_ranks$district[which(all_hrDistrict_ranks$fullDistrict_target)]
      cat(paste0("Fully-targeted districts: ", paste(fullDistrict_target_names, collapse = ", "),"\n"))

      ## calculate number of vaccines remaining for district receiving partial allocation
      partialDistrict_vacc_pop_size <- global_vacc_avail - sum(all_hrDistrict_ranks$vacc_used, na.rm = TRUE)
      ## identify highest ranked district that did not receive vaccine as the district that will receive a partial allocation
      ## NOTE: adding the criteria that pop needs to be greater than 0 to be targeted. This can occur with small polygons

      ## control flow for partial district targets
      if(sum(all_hrDistrict_ranks$fullDistrict_target, na.rm=TRUE) < nrow(all_hrDistrict_ranks)){

        partialDistrict_target_id <- all_hrDistrict_ranks$id[min(which(all_hrDistrict_ranks$vacc_used==0 & all_hrDistrict_ranks$pop>0))]
        partialDistrict_target_name <- all_hrDistrict_ranks$district[partialDistrict_target_id]
        print(partialDistrict_vacc_pop_size)
        print(partialDistrict_target_id)
        cat(sprintf("Partially-targeted district: %s (%.1f vaccinated) \n ",
                    partialDistrict_target_name, partialDistrict_vacc_pop_size))

        ## import cell-level rasters for the partial district country
        sublevel <- 2
        cntry_code <- all_hrDistrict_ranks$cntry_code[partialDistrict_target_id]
        partialDistrict_target_poly_id <- all_hrDistrict_ranks$index[partialDistrict_target_id] # index is the cntry-specific polygon ID
        cntry_code <- fix_country_name(cntry_code)
        lambda_per_cell <- R6raster$new(paste0(out_wd,cntry_code,"_lambda_per_cell.tif"))
        pop_per_cell <- R6raster$new(paste0(out_wd, cntry_code,"_",year, "_pop_per_cell.tif"))
        cases_per_cell <- R6raster$new(lambda_per_cell$raster * pop_per_cell$raster)
        cntry_shp <- get_country_sublevels(cntry_code, sublevel)

        ## get rasters for partially targeted district
        partialDistrict_lambda_rast <- extract_country_from_raster(lambda_per_cell, cntry_shp$geometry[partialDistrict_target_poly_id], partial_cover = FALSE, method = "raster")
        partialDistrict_pop_rast <- extract_country_from_raster(pop_per_cell, cntry_shp$geometry[partialDistrict_target_poly_id], partial_cover = partial_cover, method = partial_cover_method, trim_small = 0.05)
        partialDistrict_cases_rast <- R6raster$new(partialDistrict_lambda_rast$raster * partialDistrict_pop_rast$raster)

        ## allocate vaccine to the grid cells with greatest number of cases in partialDistrict at a given coverage level until there is no more
        partialDistrict_cases_summ_rast <- aggregate_raster_xlayers(partialDistrict_cases_rast, median, na.rm = TRUE)
       
        partialDistrict_lambda_summ_rast <- aggregate_raster_xlayers(partialDistrict_lambda_rast, median, na.rm = TRUE)
        ## some cells are NA, these ranks need to be last
        cell_rank <- data.frame(
            cntry_code = cntry_code,
            district = all_hrDistrict_ranks$district[partialDistrict_target_id],
            year = year,
            index = order(values(partialDistrict_lambda_summ_rast$raster), decreasing=T, na.last=TRUE),
            incidence = values(partialDistrict_lambda_summ_rast$raster)[order(values(partialDistrict_lambda_summ_rast$raster), decreasing=T, na.last=TRUE)],
            cases = values(partialDistrict_cases_summ_rast$raster)[order(values(partialDistrict_lambda_summ_rast$raster), decreasing=T, na.last=TRUE)],
            pop = values(partialDistrict_pop_rast$raster)[order(values(partialDistrict_lambda_summ_rast$raster), decreasing=T, na.last=TRUE)])

        ## get all but the last cell to target here
        cell_rank <- cell_rank %>%
        dplyr::mutate(vacc_pop_coverage = pop * as.numeric(target_coverage_prop) * doses_per_person,
                      vacc_pop_cumsum = cumsum(vacc_pop_coverage),
                      vacc_used = ifelse(vacc_pop_cumsum <= partialDistrict_vacc_pop_size, vacc_pop_coverage, 0),
                      fullCell_target = ifelse(vacc_used > 0, TRUE, FALSE)) %>%
        dplyr::mutate(fvp = vacc_used/doses_per_person,
                      prop_vaccinated = ifelse(pop==0, 0, fvp/pop))

        partialDistrict_target_cell_ids <- cell_rank[which(cell_rank$fullCell_target),]$index

        ## allocate last vaccines to highest-ranked untargeted cell, partially
        partialCell_vacc_pop_size <- partialDistrict_vacc_pop_size - sum(cell_rank$vacc_used, na.rm = TRUE)

        cat(paste0(" Full district vaccination: ", sum(all_hrDistrict_ranks$vacc_used, na.rm = TRUE),"\n"))
        cat(paste0(" Partial district vaccination: ", partialDistrict_vacc_pop_size,"\n"))
        cat(paste0(" Total vaccination: ", sum(all_hrDistrict_ranks$vacc_used, na.rm = TRUE) +
                                          partialDistrict_vacc_pop_size,"\n"))

        ## control flow for partial cell targeting
        if ((sum(cell_rank$fullCell_target, na.rm=TRUE)) < nrow(cell_rank %>% dplyr::filter(!is.na(pop)))){

          partialCell_target_cell_id <-cell_rank$index[min(which(cell_rank$vacc_used == 0 & cell_rank$pop>0), na.rm=T)]

          if(!is.na(partialCell_target_cell_id)){
            cell_rank[which(cell_rank$index==partialCell_target_cell_id),]$vacc_used <- partialCell_vacc_pop_size
            partialCell_target_prop <- (partialCell_vacc_pop_size/doses_per_person)/cell_rank[which(cell_rank$index==partialCell_target_cell_id),]$pop

            print(paste("partial cell target ID", partialCell_target_cell_id))
            print(paste("vacc pop size", partialCell_vacc_pop_size))
            print(paste("target prop", partialCell_target_prop))
          } else{
            print(paste(partialCell_vacc_pop_size, "partial cell vaccines not allocated."))
          }
          
          ## recalculate fvp and prop_vaccinated if there is partial cell targeting
          cell_rank <- cell_rank %>%
            dplyr::mutate(fvp = vacc_used/doses_per_person,
                          prop_vaccinated = ifelse(pop==0, 0, fvp/pop))
        
        } ## pc targeting

        rm(lambda_per_cell, cases_per_cell, pop_per_cell)
        gc()
      
      } ## pd targeting

      if(!dir.exists(gen_out_wd)){
        dir.create(gen_out_wd)
      }

      ## Need these outputs to create vac_template_rast by country in a separate function
      skipped_districts <- full_join(all_hrDistrict_ranks_orig, countries_ranked_by_hrDistrict, by = c("cntry_code")) %>%
        dplyr::filter(hr_district) %>%
        dplyr::mutate(skip = ifelse(district %in% skip_ls | is.na(pop), TRUE, FALSE)) %>%
        dplyr::filter(skip) %>%
        dplyr::mutate(id = NA,
                      vacc_pop_coverage = NA,
                      vacc_pop_cumsum = NA,
                      vacc_used = NA,
                      fvp = NA,
                      fullDistrict_target = NA)
      all_hrDistrict_ranks2 <- bind_rows(all_hrDistrict_ranks, skipped_districts)
      write_csv(all_hrDistrict_ranks2, paste0(gen_out_wd, scenario, "_", year, "_districtTargets.csv"))
      
      if(exists("cell_rank")){
        write_csv(cell_rank, paste0(gen_out_wd, scenario, "_", year, "_cellTargets.csv"))
        } else{
        cell_rank <- data.frame()
        }
      
      return(list(district_targets = all_hrDistrict_ranks2, cell_targets = cell_rank))

    }  else { ## no hrTargets were found
      print("No high-risk incid targets are available")
      return(list(district_targets = data.frame(), cell_targets = data.frame()))
    }
    
}


#' @name get_hrDistricts2Incid_by_country
#' @title get_hrDistricts2Incid_by_country
#' @description Rank all high risk districts in a country by incidence, where high risk districts are ones where the incidence rate is greater than 1/1000 for at least 100,000 residents or 10% of the district's population
#' @param country_year
#' @param hrThreshold1
#' @param hrThreshold2
#' @param hrThreshold3
#' @param hrThreshold4
#' @param partial_cover
#' @param partial_cover_method
#' @param out_wd
#' @return dataframe with ordered cells and incidence by country
get_hrDistricts2Incid_by_country <- function(country_year,
                                  hrThreshold1 = 1/1000,
                                  hrThreshold2 = 1/5000,
                                  hrThreshold3 = 1/10000,
                                  hrThreshold4 = 1/100000,
                                  partial_cover = FALSE,
                                  partial_cover_method = "raster",
                                  out_wd= "country_outputs/"
                                  ){

    country <- unlist(strsplit(country_year, split="--"))[1]
    year <- unlist(strsplit(country_year, split="--"))[2]
    cat(sprintf("ranking %s cells in %s \n", country, year))

    cntry <- fix_country_name(country)

    ## import pre-saved data
    lambda_per_cell <- R6raster$new(paste0(out_wd,cntry,"_lambda_per_cell.tif"))
    pop_per_district <- readRDS(paste0(out_wd,cntry,"_",year,"_pop_per_district.rds"))
    cases_per_district_samps <- readRDS(paste0(out_wd,cntry,"_",year,"_cases_per_district.rds"))
    
      print(colnames(cases_per_district_samps)[duplicated(colnames(cases_per_district_samps))])

    ## work around for duplicate ISO A2 L2 names
    cases_per_district_samps <- merge_duplicate_district_data(cases_per_district_samps)
    pop_per_district <- merge_duplicate_district_data(pop_per_district)

    inc_per_district <- apply(cases_per_district_samps/pop_per_district, 2, median)
    print("presaved data imported")
    who_region <- lookup_WorldPop_region(cntry)
    
    ## identify high risk districts
    ## summarize lambdas across samples
    cases_per_district_summ <- apply(cases_per_district_samps, 2, median, na.rm=TRUE)
    print("cases per district median")
    
    ## import existing high-risk population files if available
    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop_per_cell.tif"))){
      lambda_summ_rast <- aggregate_raster_xlayers(lambda_per_cell, median, na.rm = TRUE)
    print("lambda per cell, aggregated")
    ## import new pop version where partial_cover = FALSE
      pop_per_cell_npc <- load_population_estimates(region = who_region, lambda_per_cell, as.numeric(year))
    print("pop per cell npc, loaded")
      ## filter cells with high incidence (>1/1000)
      hr_lambda_summ_rast <- R6raster$new(raster::calc(lambda_summ_rast$raster, fun=function(x){ifelse(x > hrThreshold1, x, NA)}))
    print("hr lambda summ threshold calculated")
      ## get population for high risk cells (grab original data from non-partial cover version)
      hr_pop_per_cell <- R6raster$new(raster::mask(pop_per_cell_npc$raster, hr_lambda_summ_rast$raster, maskvalue=NA))
    print("hr pop per cell masked")

      writeRaster(hr_pop_per_cell$raster, paste0(out_wd,cntry,"_",year,"_hrpop_per_cell.tif"))
      
      rm(pop_per_cell_npc, hr_lambda_summ_rast)
      gc()

    } else {
      hr_pop_per_cell <- R6raster$new(paste0(out_wd,cntry,"_",year,"_hrpop_per_cell.tif"))
      print("hr pop per cell imported")
    }
   
    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop_per_district.rds"))){
      ## summarize hr pop by district
      hr_pop_per_district <- apply_to_all_sublevels(cntry, 
                                            hr_pop_per_cell, 
                                            ISO_level = 2, 
                                            partial_cover = partial_cover,
                                            method = partial_cover_method,
                                           trim_small = 0.05,
                                            fun=function(x){if(!all(is.na(x))){
                                              rc <- sum(x, na.rm=TRUE)
                                            } else {
                                              rc <- NA
                                            }})
      saveRDS(hr_pop_per_district, paste0(out_wd,cntry,"_",year,"_hrpop_per_district.rds"))
    print("hr pop per district saved")
    } else{
      hr_pop_per_district <- readRDS(paste0(out_wd,cntry,"_",year,"_hrpop_per_district.rds"))
    print("hr pop per district imported")
    }

  ## work around for duplicate ISO A2 L2 names
  hr_pop_per_district <- merge_duplicate_district_data(hr_pop_per_district)

    ## import second-level high-risk population files if available
    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop5k_per_cell.tif"))){
      lambda_summ_rast <- aggregate_raster_xlayers(lambda_per_cell, median, na.rm = TRUE)
    print("lambda per cell, aggregated")
    ## import new pop version where partial_cover = FALSE
      pop_per_cell_npc <- load_population_estimates(region = who_region, lambda_per_cell, as.numeric(year))
    print("pop per cell npc, loaded")
      ## filter cells with high incidence (>1/1000)
      hr2_lambda_summ_rast <- R6raster$new(raster::calc(lambda_summ_rast$raster, fun=function(x){ifelse(x > hrThreshold2, x, NA)}))
    print("hr2 lambda summ threshold calculated")
      ## get population for high risk cells (grab original data from non-partial cover version)
      hr2_pop_per_cell <- R6raster$new(raster::mask(pop_per_cell_npc$raster, hr2_lambda_summ_rast$raster, maskvalue=NA))
    print("hr2 pop per cell masked")

      writeRaster(hr2_pop_per_cell$raster, paste0(out_wd,cntry,"_",year,"_hrpop5k_per_cell.tif"))
      
      rm(pop_per_cell_npc, hr2_lambda_summ_rast)
      gc()
    } else{
      hr2_pop_per_cell <- R6raster$new(paste0(out_wd,cntry,"_",year,"_hrpop5k_per_cell.tif"))
    print("hr2 pop per cell imported")
    }

    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop5k_per_district.rds"))){
      ## summarize hr pop by district
      hr2_pop_per_district <- apply_to_all_sublevels(cntry, 
                                            hr2_pop_per_cell, 
                                            ISO_level = 2, 
                                            partial_cover = partial_cover,
                                            method = partial_cover_method,
                                           trim_small = 0.05,
                                            fun=function(x){if(!all(is.na(x))){
                                              rc <- sum(x, na.rm=TRUE)
                                            } else {
                                              rc <- NA
                                            }})
      saveRDS(hr2_pop_per_district, paste0(out_wd,cntry,"_",year,"_hrpop5k_per_district.rds"))
    print("hr2 pop per district saved")
    } else{
      hr2_pop_per_district <- readRDS(paste0(out_wd,cntry,"_",year,"_hrpop5k_per_district.rds"))
    print("hr2 pop per district imported")
    }

    ## work around for duplicate ISO A2 L2 names
    hr2_pop_per_district <- merge_duplicate_district_data(hr2_pop_per_district)

    ## import third-level high-risk population files if available
    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop10k_per_cell.tif"))){
      lambda_summ_rast <- aggregate_raster_xlayers(lambda_per_cell, median, na.rm = TRUE)
    print("lambda per cell, aggregated")
    ## import new pop version where partial_cover = FALSE
      pop_per_cell_npc <- load_population_estimates(region = who_region, lambda_per_cell, as.numeric(year))
    print("pop per cell npc, loaded")
      ## filter cells with high incidence (>1/10000)
      hr3_lambda_summ_rast <- R6raster$new(raster::calc(lambda_summ_rast$raster, fun=function(x){ifelse(x > hrThreshold3, x, NA)}))
    print("hr3 lambda summ threshold calculated")
      ## get population for high risk cells (grab original data from non-partial cover version)
      hr3_pop_per_cell <- R6raster$new(raster::mask(pop_per_cell_npc$raster, hr3_lambda_summ_rast$raster, maskvalue=NA))
    print("hr3 pop per cell masked")

      writeRaster(hr3_pop_per_cell$raster, paste0(out_wd,cntry,"_",year,"_hrpop10k_per_cell.tif"))
      
      rm(pop_per_cell_npc, hr3_lambda_summ_rast)
      gc()
    } else{
      hr3_pop_per_cell <- R6raster$new(paste0(out_wd,cntry,"_",year,"_hrpop10k_per_cell.tif"))
  print("hr3 pop per cell imported")
    }

    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop10k_per_district.rds"))){
      ## summarize hr pop by district
      hr3_pop_per_district <- apply_to_all_sublevels(cntry, 
                                            hr3_pop_per_cell, 
                                            ISO_level = 2, 
                                            partial_cover = partial_cover,
                                            method = partial_cover_method,
                                           trim_small = 0.05,
                                            fun=function(x){if(!all(is.na(x))){
                                              rc <- sum(x, na.rm=TRUE)
                                            } else {
                                              rc <- NA
                                            }})
      saveRDS(hr3_pop_per_district, paste0(out_wd,cntry,"_",year,"_hrpop10k_per_district.rds"))
  print("hr3 pop per district saved")
    } else{
      hr3_pop_per_district <- readRDS(paste0(out_wd,cntry,"_",year,"_hrpop10k_per_district.rds"))
  print("hr3 pop per district imported")
    }

    ## work around for duplicate ISO A2 L2 names
    hr3_pop_per_district <- merge_duplicate_district_data(hr3_pop_per_district)

    ## import fourth-level high-risk population files if available
    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop100k_per_cell.tif"))){
      lambda_summ_rast <- aggregate_raster_xlayers(lambda_per_cell, median, na.rm = TRUE)
    print("lambda per cell, aggregated")
    ## import new pop version where partial_cover = FALSE
      pop_per_cell_npc <- load_population_estimates(region = who_region, lambda_per_cell, as.numeric(year))
    print("pop per cell npc, loaded")
      ## filter cells with high incidence (>1/10000)
      hr4_lambda_summ_rast <- R6raster$new(raster::calc(lambda_summ_rast$raster, fun=function(x){ifelse(x > hrThreshold4, x, NA)}))
    print("hr4 lambda summ threshold calculated")
      ## get population for high risk cells (grab original data from non-partial cover version)
      hr4_pop_per_cell <- R6raster$new(raster::mask(pop_per_cell_npc$raster, hr4_lambda_summ_rast$raster, maskvalue=NA))
    print("hr4 pop per cell masked")

      writeRaster(hr4_pop_per_cell$raster, paste0(out_wd,cntry,"_",year,"_hrpop100k_per_cell.tif"))
      
      rm(pop_per_cell_npc, hr4_lambda_summ_rast)
      gc()
    } else{
      hr4_pop_per_cell <- R6raster$new(paste0(out_wd,cntry,"_",year,"_hrpop100k_per_cell.tif"))
  print("hr4 pop per cell imported")
    }

    if(!file.exists(paste0(out_wd,cntry,"_",year,"_hrpop100k_per_district.rds"))){
      ## summarize hr pop by district
      hr4_pop_per_district <- apply_to_all_sublevels(cntry, 
                                            hr4_pop_per_cell, 
                                            ISO_level = 2, 
                                            partial_cover = partial_cover,
                                            method = partial_cover_method,
                                           trim_small = 0.05,
                                            fun=function(x){if(!all(is.na(x))){
                                              rc <- sum(x, na.rm=TRUE)
                                            } else {
                                              rc <- NA
                                            }})
      saveRDS(hr4_pop_per_district, paste0(out_wd,cntry,"_",year,"_hrpop100k_per_district.rds"))
  print("hr4 pop per district saved")
    } else{
      hr4_pop_per_district <- readRDS(paste0(out_wd,cntry,"_",year,"_hrpop100k_per_district.rds"))
  print("hr4 pop per district imported")
    }

    ## work around for duplicate ISO A2 L2 names
    hr4_pop_per_district <- merge_duplicate_district_data(hr4_pop_per_district)
    
    ## highest risk districts
    hr1_df <- tbl_df(data.frame(cntry_code = cntry, 
                                  year = year, 
                                  district = colnames(pop_per_district)[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  pop = pop_per_district[1, order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  hr_pop = hr_pop_per_district[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  cases = cases_per_district_summ[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  incidence = inc_per_district[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  index = order(inc_per_district, decreasing=T, na.last=TRUE))) %>%
      dplyr::mutate(hr_prop = ifelse(pop > 0, hr_pop/pop, 0),
                    hr_district = ifelse(hr_prop > 0.1 | hr_pop > 100000, TRUE, FALSE),
                    hr_level = 1)
    hr_dnames <- hr1_df %>% dplyr::filter(!is.na(hr_pop)) %>% dplyr::distinct(district) %>% unlist

    ## second-level high risk districts
    hr2_df <- tbl_df(data.frame(cntry_code = cntry, 
                                  year = year, 
                                  district = colnames(pop_per_district)[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  pop = pop_per_district[1, order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  hr_pop = hr2_pop_per_district[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  cases = cases_per_district_summ[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  incidence = inc_per_district[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  index = order(inc_per_district, decreasing=T, na.last=TRUE))) %>%
      dplyr::mutate(hr_prop = ifelse(pop > 0, hr_pop/pop, 0),
                    hr_district = ifelse(hr_prop > 0.1 | hr_pop > 100000, TRUE, FALSE),
                    hr_level = 2) %>%
      dplyr::filter(!(district %in% hr_dnames)) ## rm hrDistrict level 1
    hr_dnames2 <- hr2_df %>% dplyr::filter(!is.na(hr_pop)) %>% dplyr::distinct(district) %>% unlist


    ## third-level high risk districts
    hr3_df <- tbl_df(data.frame(cntry_code = cntry, 
                                  year = year, 
                                  district = colnames(pop_per_district)[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  pop = pop_per_district[1, order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  hr_pop = hr3_pop_per_district[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  cases = cases_per_district_summ[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  incidence = inc_per_district[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  index = order(inc_per_district, decreasing=T, na.last=TRUE))) %>%
      dplyr::mutate(hr_prop = ifelse(pop > 0, hr_pop/pop, 0),
                    hr_district = ifelse(hr_prop > 0.1 | hr_pop > 100000, TRUE, FALSE),
                    hr_level = 3) %>%
      dplyr::filter(!(district %in% hr_dnames) & !(district %in% hr_dnames2)) ## rm hrDistrict level 1 & 2
    hr_dnames3 <- hr3_df %>% dplyr::filter(!is.na(hr_pop)) %>% dplyr::distinct(district) %>% unlist

    ## fourth-level high risk districts
    hr4_df <- tbl_df(data.frame(cntry_code = cntry, 
                                  year = year, 
                                  district = colnames(pop_per_district)[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  pop = pop_per_district[1, order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  hr_pop = hr4_pop_per_district[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  cases = cases_per_district_summ[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  incidence = inc_per_district[order(inc_per_district, decreasing=T, na.last=TRUE)], 
                                  index = order(inc_per_district, decreasing=T, na.last=TRUE))) %>%
      dplyr::mutate(hr_prop = ifelse(pop > 0, hr_pop/pop, 0),
                    hr_district = ifelse(hr_prop > 0.1 | hr_pop > 100000, TRUE, FALSE),
                    hr_level = 4) %>%
      dplyr::filter(!(district %in% hr_dnames) & !(district %in% hr_dnames2) & !(district %in% hr_dnames3)) ## rm hrDistrict level 1 & 2 & 3

    hr_df <- bind_rows(hr1_df, hr2_df, hr3_df, hr4_df)
    
  rm(lambda_per_cell, hr_pop_per_cell)
  gc()

  return(hr_df)
}


#' @name skip_district_targets
#' @title skip_district_targets
#' @description Skip districts that were targets in specific years (to simulate the rotation of mass vaccination campaigns in different years)
#' @param skipYears
#' @param scenario
#' @param out_wd
#' @return list of fully targeted districts from skipYears
skip_district_targets <- function(skipYears, 
                                  scenario,
                                  out_wd = "country_targets/"){
  
    skip_targets_files <- paste0(out_wd, scenario, "_", skipYears, "_districtTargets.csv")

    skip_ls <- lapply(skip_targets_files, function(x){
        if(file.exists(x)) {
          indat <- read_csv(x) 
          if("fullDistrict_target" %in% names(indat)){
            skips <- indat  %>% dplyr::filter(fullDistrict_target) %>%
            dplyr::select(district) %>%
            unlist %>% unname
          } else{
            print(paste(x, "does not have full district targets"))
            skips <- c(NA)
          }
        } else{
          print(paste(x, "does not exist"))
          skips <- c(NA)
        }
        return(skips)
      }) %>% unlist %>% unique 
    skip_ls <- skip_ls[!is.na(skip_ls)]

    return(skip_ls)
}


#' @name skip_hrDistrict_targets2
#' @title skip_hrDistrict_targets2
#' @description Skip high incidence districts targeted in specific years (to simulate the rotation of mass vaccination campaigns in different years), following country rank, specifically for the rate-logistics and case-logistics models
#' @param skipYears
#' @param scenario
#' @param out_wd
#' @return list of fully targeted high incidence districts from skipyears
skip_hrDistrict_targets2 <- function(skipYears, scenario, out_wd = "country_targets/"){
  
  skip_targets_files <- paste0(out_wd, scenario, "_", skipYears, "_districtTargets.csv")
  skip_ls <- lapply(skip_targets_files, function(x){
    if(file.exists(x)){
      skip_summ <- read_csv(x) %>% 
        dplyr::filter(fullDistrict_target) %>%
        dplyr::select(district) %>%
        unlist 
    } else{
      print(paste(x, "does not exist"))
      skip_summ <- c(NA)
    }
    return(skip_summ)
  }) %>% unlist %>% unique 
  skip_ls <- skip_ls[!is.na(skip_ls)]

  return(skip_ls)
}


#' @name merge_duplicate_district_data
#' @title merge_duplicate_district_data
#' @description Workaround for duplicate ISO A2 L2 names, merge duplicate districts (e.g., pop_per_district)
#' @param district_object
#' @return dataframe with merged district data
merge_duplicate_district_data <- function(district_object){
  dupnames <- colnames(district_object)[duplicated(colnames(district_object))]
  if(length(dupnames)>0){
    new_do <- district_object
    for(i in 1:length(dupnames)){
      nm <- dupnames[i]
      ix <- which(colnames(district_object)==dupnames[i])
      new_col <- apply(district_object[,ix], 1, sum, na.rm=TRUE)
      new_do <- cbind(district_object[,-ix], new_col)
      colnames(new_do) <- c(colnames(new_do)[-length(colnames(new_do))], nm)
    }
    district_object <- new_do
  }
  return(district_object)
}


#' @name identify_optimumWashPullan_targets
#' @title identify_optimumWashPullan_targets
#' @description Allocate vaccines according to districts with the least access to improved water and sanitation, independent of country
#' @param wash_district_df 
#' @param cntryCode_ls
#' @param scenario
#' @param vacc_inputs
#' @param doses_per_person
#' @param rotate_every_nYears
#' @param gen_out_wd
#' @return list with one dataframe, one for district-level (fullDistrict & partialDistrict) targets 
identify_optimumWashPullan_targets <- function(wash_district_df,
                              cntryCode_ls,
                              scenario,
                              vacc_inputs,
                              method_coverage_prop,
                              doses_per_person = 2,
                              rotation = 3,
                              gen_in_wd = "country_outputs/",
                              gen_out_wd = "country_targets/"
                              ){

    cat(sprintf("allocate vaccines with %s targeting by cases for %s: %s \n", scenario, "allyrs", paste(cntryCode_ls, collapse=", ")))

    covEstimate <- get_vacc_coverage() %>% dplyr::filter(method == method_coverage_prop) %>% dplyr::select(estimate) %>% unlist
    vacc_inputs <- vacc_inputs %>% dplyr::mutate(coverage = covEstimate)

    all_targets <- map_dfr(1:nrow(vacc_inputs), function(i){
      colname <- get_washname(scenario)
      target_coverage_prop <- as.numeric(vacc_inputs$coverage[i])
      global_vacc_avail <- as.numeric(vacc_inputs$vaccines[i])
      yr = vacc_inputs$year[i]

      all_district_ranks <- rank_washPullan_by_alldistricts(countryCode_ls = cntryCode_ls, scenario = scenario, wash_df = wash_district_df, year = yr, rotate_every_nYears = rotation, out_wd = gen_out_wd, in_wd = gen_in_wd) %>%
        dplyr::filter(!skip) %>%
        dplyr::arrange(desc(!!sym(colname))) %>%
        dplyr::mutate(id = seq_along(pop),
                      vacc_pop_coverage = pop * target_coverage_prop * doses_per_person,
                      vacc_pop_cumsum = cumsum(vacc_pop_coverage),
                      vacc_used = ifelse(vacc_pop_cumsum <= global_vacc_avail, vacc_pop_coverage, 0),
                      fvp = vacc_used/doses_per_person,
                      fullDistrict_target = ifelse(vacc_used > 0, TRUE, FALSE))

      print("alldistricts")
      print(head(all_district_ranks))
      
      fullDistrict_target_names <- all_district_ranks$district[which(all_district_ranks$fullDistrict_target)]
      cat(paste0("Fully-targeted districts: ", paste(fullDistrict_target_names, collapse = ", "),"\n"))

      ## calculate number of vaccines remaining for district receiving partial allocation
      partialDistrict_vacc_pop_size <- global_vacc_avail - sum(all_district_ranks$vacc_used, na.rm = TRUE)
      
      ## identify highest ranked district that did not receive vaccine as the district that will receive a partial allocation
      ## NOTE: adding the criteria that pop needs to be greater than 0 to be targeted. This can occur with small polygons
      partialDistrict_target_id <- all_district_ranks$id[min(which(all_district_ranks$vacc_used==0 & all_district_ranks$pop>0))]
      partialDistrict_target_name <- all_district_ranks$district[partialDistrict_target_id]
      print(partialDistrict_target_id)
      cat(sprintf("Partially-targeted district: %s (%.1f vaccinated) \n ",
                  partialDistrict_target_name, partialDistrict_vacc_pop_size))

      ## insert partial vacc district values into orig data frame
      all_district_ranks[which(all_district_ranks$id == partialDistrict_target_id),]$vacc_used <- partialDistrict_vacc_pop_size
      all_district_ranks[which(all_district_ranks$id == partialDistrict_target_id),]$fvp <- partialDistrict_vacc_pop_size/doses_per_person

      ## Need these outputs to create vac_template_rast by country in a separate function
      write_csv(all_district_ranks, paste0(gen_out_wd, "/", scenario, "_", yr, "_districtTargets.csv"))

    })
    
    return(list(district_targets = all_targets))
}


#' @name rank_washPullan_by_alldistricts
#' @title rank_washPullan_by_alldistricts
#' @description Rank all districts in Africa by water, sanitation, or joint coverage according to Pullan et al. 2012
#' @param countryCode_ls
#' @param scenario
#' @param wash_df
#' @param year
#' @param rotate_every_nYears
#' @param in_wd
#' @param out_wd
#' @return dataframe with ordered wash coverage by country
rank_washPullan_by_alldistricts <- function(countryCode_ls,
                                  scenario,
                                  wash_df, 
                                  year,
                                  rotate_every_nYears = 3,
                                  in_wd = "country_outputs/",
                                  out_wd = "country_targets/"
                                  ){

  ## which WASH measure is used?
  colname <- get_washname(scenario)

  ## get list of fully targeted districts in previous x years
  skipYears <- year-(1:(rotate_every_nYears-1))
  skip_ls <- skip_district_targets(skipYears, scenario, out_wd)

  ## grab population in appropriate model year
  pops_per_district <- map_dfr(countryCode_ls, function(cntry){
    inpop <- readRDS(file = paste0(in_wd,cntry,"_",year,"_pop_per_district.rds"))
    pop <- data.frame(district = colnames(inpop), pop = inpop[1,], row.names = NULL, stringsAsFactors = FALSE)
    return(pop)
  })
  
  print(head(pops_per_district))
  print(head(wash_df))
  wash_df2 <- wash_df %>%
    dplyr::select(-pop) %>%
    dplyr::mutate(year = year) %>%
    dplyr::filter(cntry_code %in% countryCode_ls) %>%
    dplyr::select(cntry_code, year, district, wat, san, wash) %>%
    right_join(pops_per_district, by = c("district")) %>%
    group_by(cntry_code) %>%
    dplyr::arrange(desc(!!sym(colname))) %>%
    dplyr::mutate(index = seq_along(pop)) %>%
    ungroup %>%
    dplyr::mutate(skip = ifelse(district %in% skip_ls | is.na(pop), TRUE, FALSE))
  print(head(wash_df2))

  return(wash_df2)
}


#' @name run_country_scenario
#' @title run_country_scenario
#' @description Run model by country
#' @param country
#' @param scenario
#' @param years year of study projection
#' @param alloc_strategy
#' @param target_inc_threshold
#' @param mu population turnover rate
#' @param ve_direct vaccine efficacy function
#' @param indirect_mult indirect vaccine protection multiplier
#' @param secular_trend_multiplier generate_flatline_multiplier() projects constant risk in future years 
#' @param track_vac
#' @param use_partial_cover
#' @param method
#' @param in_wd
#' @param out_wd
#' @param cf_wd
#' @param targets_wd
#' @param gen_out_wd
#' @return dataframe with expected cases and vaccinated population per year
run_country_scenario <- function(country,
                                 scenario,
                                 years = 2018:2030,  
                                 alloc_strategy,                 
                                 target_inc_threshold = 1/1000,
                                 mu = 1/65,
                                 ve_direct = generate_pct_protect_function(),
                                 indirect_mult = generate_indirect_incidence_mult(),
                                 secular_trend_multiplier = generate_flatline_multiplier(),
                                 track_vac = TRUE,
                                 use_partial_cover = FALSE,
                                 method = "raster",
                                 in_wd="data/",
                                 out_wd="country_outputs/",
                                 cf_wd = "generated_outputs/cf/",
                                 targets_wd="country_targets/",
                                 gen_out_wd="generated_outputs/"
                                 ){
    
    cntry <- fix_country_name(country)

    ## import lambda raster for force of infection input and to use as a dummy rasterStack
    lambda_fn <- paste0(out_wd,cntry,"_lambda_per_cell.tif")
    lambda_per_cell <- R6raster$new(lambda_fn)

    ## Initiate rasterLayer to begin target rasterStack
    dummy_rasterLayer <- R6raster$new(paste0(out_wd, cntry, "_", years[1], "_pop_per_cell.tif"))
    vac_rast <- R6raster$new(dummy_rasterLayer$raster[[1]])

    ## loop over all years of interest
    for (i in 1:length(years)){
      
      tmp_vac_rast <- create_vacc_allocation_raster(scenario = scenario,
                                                    alloc_strategy = alloc_strategy, 
                                                    country = country,
                                                    gen_out_wd = targets_wd,
                                                    out_wd = out_wd,
                                                    year = years[i]) 
      
      vac_rast <- R6raster$new(raster::stack(vac_rast$raster, tmp_vac_rast$raster))
      ## clear up memory by removing unneeded rasters
      rm(tmp_vac_rast)
      gc()

    }
    
    ## clear up memory by removing unneeded rasters
    rm(dummy_rasterLayer)
    gc()
    
    ## Delete the first raster in the stack because it was dummy data
    vac_rast <- R6raster$new(vac_rast$raster[[-1]])

    pop_rast <- R6raster$new(sapply(years,
                       function(x) R6raster$new(paste0(out_wd, cntry, "_", x, "_pop_per_cell.tif"))$raster) %>%
    stack)

    ## create dummy rasterStack
    sus_rast <- R6raster$new(vac_rast$raster)
    track_prop_immune <- data.frame()

    ## skip the modeling part if all of the outputs exist
    if (!file.exists(paste0(cf_wd, cntry, "_cf_cases_samples.csv")) | 
      !file.exists(paste0(gen_out_wd, cntry, "_", scenario, "_cases_samples.csv")) | !file.exists(paste0(gen_out_wd, cntry, "_", "sus_", scenario, ".tif"))){
        
      print(paste("inverse mu", 1/mu)) 
      ## now step through each layer and update sus_rast with the number
      ## of people susceptible each year
      for (t in 1:length(years)){ #by year
          ## create a template for multiplying
          tmp <- rep(1, length(values(sus_rast$raster[[1]])))

          for (k in 1:t){
              pkt <- values(pop_rast$raster[[k]])*(1-(t-k)*mu)/values(pop_rast$raster[[t]])
              prob_still_protected <- values(vac_rast$raster[[k]]) * c(ve_direct(t-k+1)) * pkt
              tmp <- tmp * (1 - prob_still_protected)
              
              ## track proportion of full population immune for cases averted calculation
              prop_protected <- sum(prob_still_protected * values(pop_rast$raster[[t]]), na.rm=T)/sum(values(pop_rast$raster[[t]]), na.rm=T)
              track_prop_immune <- bind_rows(track_prop_immune,
                                             c(year=years[t], campaignYear=years[k], prop_immune_in_year=prop_protected))
          }

          values(sus_rast$raster[[t]]) <- tmp
      }

      writeRaster(sus_rast$raster, paste0(gen_out_wd, cntry, "_", "sus_", scenario, ".tif"))
      
      ## Save case and counterfactual samples by year so cumulative health impacts can be calculated with error bars:
      nlayers <- dim(lambda_per_cell$raster)[3]
      cases_samples_vac <- data.frame(matrix(, ncol=length(years), nrow=nlayers))
      cases_samples_noVac <-  data.frame(matrix(, ncol=length(years), nrow=nlayers))
      
      ## The for loop will be slower for processing but more memory efficient than lapply
      for(t in 1:length(years)) {
        ## if there are no targets in any years, skip the vacc scenario and save the no vacc scenario twice
        if(sum(values(vac_rast$raster), na.rm = TRUE) > 0) {
          
          print(paste("Starting vacc scenario", cntry, years[t]))
          if(file.exists(paste0(gen_out_wd, cntry, "_", years[t], "_", scenario, ".tif"))){
            ec_withVac_rast <- R6raster$new(paste0(gen_out_wd, cntry, "_", years[t], "_", scenario, ".tif"))

          } else {
            indirect_tmp <- R6raster$new(sus_rast$raster[[t]])
            values(indirect_tmp$raster) <- indirect_mult(1-values(sus_rast$raster[[t]]))

            ec_withVac_rast <- R6raster$new(sus_rast$raster[[t]]*
              pop_rast$raster[[t]]*
              lambda_per_cell$raster*
              indirect_tmp$raster*
              secular_trend_multiplier(year = years[t]))

            print(paste("Saving vacc scenario", cntry, years[t]))
            writeRaster(ec_withVac_rast$raster, paste0(gen_out_wd, cntry, "_", years[t], "_", scenario, ".tif"))

            rm(indirect_tmp)
            gc()
          }
          
          print(paste("Processing case samples", cntry, years[t]))
          
          cases_samples_vac_yr <- aggregate_raster_xcells(ec_withVac_rast, sum, na.rm = TRUE)
          cases_samples_vac[,t] <- unname(cases_samples_vac_yr)
          names(cases_samples_vac)[t] <- paste0("x", years[t])


          ## clear up memory by removing unneeded rasters
          rm(ec_withVac_rast, cases_samples_vac_yr)
          gc()
          
          print(paste("Ending vacc scenario", cntry, years[t]))
        } else {
          print(paste("Skipping vacc scenario, no targets", cntry, years[t]))
        }
          
        print(paste("Starting no vacc scenario", cntry, years[t]))
        if(file.exists(paste0(cf_wd, cntry, "_", years[t], "_cf.tif"))) {
          ec_noVac_rast <- R6raster$new(paste0(cf_wd, cntry, "_", years[t], "_cf.tif"))
        } else{ 
          ec_noVac_rast <- R6raster$new(pop_rast$raster[[t]]*
            lambda_per_cell$raster*
            secular_trend_multiplier(year = years[t]))

          print(paste("Saving no vacc scenario", cntry, years[t]))
          writeRaster(ec_noVac_rast$raster, paste0(cf_wd, cntry, "_", years[t], "_cf.tif"))
        }

        ## if no vaccination, copy no vacc scenarios for vacc scenarios
        if(sum(values(vac_rast$raster), na.rm = TRUE) == 0) {
          print(paste("Duplicate saving no vacc scenario as vacc scenario", cntry, years[t]))
          if(!file.exists(paste0(gen_out_wd, cntry, "_", years[t], "_", scenario, ".tif"))) {
            writeRaster(ec_noVac_rast$raster, paste0(gen_out_wd, cntry, "_", years[t], "_", scenario, ".tif"))
          }
        }
          
        print(paste("Processing no vacc samples", cntry, years[t]))
        cases_samples_noVac_yr <- aggregate_raster_xcells(ec_noVac_rast, sum, na.rm = TRUE)
        cases_samples_noVac[,t] <- unname(cases_samples_noVac_yr)
        names(cases_samples_noVac)[t] <- paste0("x", years[t])

        rm(ec_noVac_rast, cases_samples_noVac_yr)
        gc()  

        print(paste("Ending no vacc scenario", cntry, years[t])) 
      } ## end for loop

      ## clear up memory by removing lambda raster
      rm(lambda_per_cell, sus_rast)
      gc()

      ## if no vaccination, copy no vacc scenarios for vacc scenarios
      if(sum(values(vac_rast$raster), na.rm = TRUE) == 0) {
        cases_samples_vac <- cases_samples_noVac
      }  

      ## export cases and cf samples to csv
      cases_samples_vac <- cases_samples_vac %>%
        dplyr::mutate(cntry_code=cntry, allocation_strategy=alloc_strategy, scen=scenario)  
      write_csv(cases_samples_vac, paste0(gen_out_wd, cntry, "_", scenario, "_cases_samples.csv"))

      if(!file.exists(paste0(cf_wd, cntry, "_cf_cases_samples.csv"))) {
        cases_samples_noVac <- cases_samples_noVac %>%
            dplyr::mutate(cntry_code=cntry, allocation_strategy=alloc_strategy, scen=scenario)
        write_csv(cases_samples_noVac, paste0(cf_wd, cntry, "_cf_cases_samples.csv"))
      }

    } else {

      print(paste("Loading case and cf samples from file", cntry, scenario))
      cases_samples_vac <- read_csv(paste0(gen_out_wd, cntry, "_", scenario, "_cases_samples.csv"))
      cases_samples_noVac <- read_csv(paste0(cf_wd, cntry, "_cf_cases_samples.csv"))
    }

    cases_samples_vac <- cases_samples_vac %>% dplyr::select(contains("x20"))
    cases_samples_noVac <- cases_samples_noVac %>% dplyr::select(contains("x20"))

    ## process country-level summary from samples
    print(paste("Processing country-level summary from samples"))
    cases_country_vac <- map_dfr(1:ncol(cases_samples_vac), function(x){
        cases_country_vac_yr <- quantile(as.numeric(unlist(cases_samples_vac[,x])), prob=c(.5, .025, .975)) %>%
          t %>%
          tidy %>%
          rename(median_cases=X50., ci_lower_cases=X2.5., ci_upper_cases=X97.5.)
        return(cases_country_vac_yr)
      }) %>%
      dplyr::mutate(cntry_code=cntry, year=years)
    rm(cases_samples_vac)
    gc()

    cases_country_noVac <- map_dfr(1:ncol(cases_samples_noVac), function(x){
      cases_country_noVac_yr <- quantile(as.numeric(unlist(cases_samples_noVac[,x])), prob=c(.5, .025, .975)) %>%
        t %>%
        tidy %>%
        rename(median_cases=X50., ci_lower_cases=X2.5., ci_upper_cases=X97.5.)
      return(cases_country_noVac_yr)
      }) %>%
      dplyr::rename(median_cases_cf = median_cases, ci_lower_cases_cf = ci_lower_cases, ci_upper_cases_cf = ci_upper_cases) %>%
      dplyr::mutate(cntry_code=cntry, year=years)
    rm(cases_samples_noVac)
    gc()

    ## track size of vaccinated population by year
    if(track_vac){
      print("Tracking vaccinated population")
      vac_track_mat <- R6raster$new(vac_rast$raster * pop_rast$raster)
      cases_country_vac <- inner_join(cases_country_vac,
                                      data.frame(year=years,
                                                 fvp=sapply(1:(dim(vac_track_mat$raster)[3]), function(x) sum(values(vac_track_mat$raster[[x]]), na.rm=TRUE))))
    }
      
    summary_cases_out <- full_join(cases_country_vac, cases_country_noVac,
                                       by = c("cntry_code", "year")) %>%
      dplyr::mutate(runDate = Sys.Date()) %>%
      dplyr::select(runDate, cntry_code, year, fvp, median_cases, median_cases_cf, ci_lower_cases, ci_upper_cases, ci_lower_cases_cf, ci_upper_cases_cf) 
    
    ## calculate health outcomes from cases averted raster
    print(paste("Saving summary health outcomes", cntry, scenario))
    summary_outcomes_out <- run_country_outputs(country = cntry, scenario = scenario, gen_out_wd = gen_out_wd, cf_wd = cf_wd)
    summary_out <- left_join(summary_cases_out, summary_outcomes_out, by = c("cntry_code", "year"))
    write_csv(summary_out, paste0(gen_out_wd, cntry, "_", scenario, "_summary_out.csv"))

    rm(vac_rast, pop_rast, vac_track_mat, summary_cases_out, summary_outcomes_out)
    gc()
    return(summary_out)

}


#' @name create_vacc_allocation_raster
#' @title create_vacc_allocation_raster
#' @description Allocate vaccines according to targets listed in the previously exported district or cell target spreadsheets
#' @param scenario
#' @param alloc_strategy
#' @param country
#' @param year
#' @param plot_me
#' @param out_wd
#' @param gen_out_wd
#' @return rasterLayer with the proportion of the population allocated with vaccines in a given country and year
create_vacc_allocation_raster <- function(scenario,
                              alloc_strategy, 
                              country,
                              year,
                              plot_me = FALSE,
                              out_wd = "country_outputs/",
                              gen_out_wd = "country_targets/"
                              ){

  sublevel <- 2 ## district-level (ISO A2 L2 administrative units)
  cntry <- fix_country_name(country)
  cntry_shp <- get_country_sublevels(country, sublevel)
  pop_per_cell <- R6raster$new(paste0(out_wd, cntry, "_", year, "_pop_per_cell.tif"))
  vac_template_rast <- R6raster$new(pop_per_cell$raster)
  ## Put the allocations into vac_template_rast as proportion of pop covered in that cell by vaccine
  ## create raster template from previously imported pop raster
  values(vac_template_rast$raster) <- ifelse(!is.na(values(vac_template_rast$raster)),0,NA)
  rm(pop_per_cell)
  gc()

  ## grab district targets for country (if they exist)
  if(file.exists(paste0(gen_out_wd, scenario, "_", year, "_districtTargets.csv")) & !grepl("allDistrict_", scenario)) {
    district_targets <- read_csv(paste0(gen_out_wd, scenario, "_", year, "_districtTargets.csv")) %>%
      dplyr::filter(cntry_code == cntry & fullDistrict_target == TRUE)
    ## calculate target coverage proportion from districtTargets file
    target_coverage_prop <- district_targets[1,]$fvp/district_targets[1,]$pop 
  } else if(file.exists(paste0(gen_out_wd, scenario, "_", year, "_districtTargets.csv")) & grepl("allDistrict_", scenario)) {
    district_targets <- read_csv(paste0(gen_out_wd, scenario, "_", year, "_districtTargets.csv")) %>%
      dplyr::filter(cntry_code == cntry)
    ## grab target coverage proportion from districtTargets file
    target_coverage_prop <- district_targets[1,]$prop_vaccinated
  } else {
    district_targets <- data.frame()
  }
    
  ## grab cell targets for the country (if they exist)
  if(file.exists(paste0(gen_out_wd, scenario, "_", year, "_cellTargets.csv"))){
    cell_targets <- read_csv(paste0(gen_out_wd, scenario, "_", year, "_cellTargets.csv")) %>%
      dplyr::filter(cntry_code == cntry & prop_vaccinated > 0)
  } else{
    cell_targets <- data.frame()
  }
  

  if (prod(dim(district_targets))>0 | prod(dim(cell_targets))>0) {
    ## proportion of population vaccinated in fully covered districts = level of vacc coverage
    if (prod(dim(district_targets))>0){
    print(sprintf("adding district targets for %s %s", country, year))
    
    district_target_poly_ids <- district_targets$index

    for (i in district_target_poly_ids){
        vac_template_rast$raster <- raster::mask(vac_template_rast$raster,
                                  as(cntry_shp$geometry[i], 'Spatial'),
                                  inverse=TRUE,
                                  updatevalue=target_coverage_prop)
      }
    }

    ## cell targeting
    if (prod(dim(cell_targets))>0){
      print(sprintf("adding cell targets for %s %s", country, year))
      cell_target_cell_ids <- cell_targets$index
      
      for (j in 1:nrow(cell_targets)){
          values(vac_template_rast$raster)[cell_targets$index[j]] <- cell_targets$prop_vaccinated[j]
      }
    }

  } else{
    print(sprintf("no targets for %s %s", country, year))
    # return(vac_template_rast) # raster should only have 0 or NA
  }

  ## check: overlay target cells on country map
  if(plot_me){
      plot(vac_template_rast$raster,main=paste(cntry, year))
      plot(cntry_shp$geometry)
      if (prod(dim(cell_targets))>0){
      plot(cntry_shp$geometry[cell_target_cell_ids], add=T, col="blue")
      }
      if (prod(dim(district_targets))>0){
      plot(cntry_shp$geometry[district_target_poly_ids], add=T, col="red")
     }
      
  }
  return(vac_template_rast)
} 


#' @name get_vacc_cost_data
#' @title get_vacc_cost_data
#' @description Import OCV campaign cost surveys spreadsheet and process delivery, procurement, and total costs per FVP for a single parameter set
#' @param data_wd
#' @param sensitivity
#' @return 
get_vacc_cost_data <- function(data_wd = "data/",
                                sensitivity = "base"){

    cost_dat <- read_xlsx(paste0(data_wd, "vacc_costs.xlsx"), sheet = 1, col_names = TRUE, col_types = "text", na = c("NA", "TBD")) %>%
      dplyr::mutate(cost_survey = as.logical(cost_survey), incl = as.logical(incl)) %>%
      dplyr::filter(who_region == "AFR" & !is.na(vacc_delivery_cost_per_FVP) & incl) %>%
      rowwise %>%
      dplyr::mutate(orig_delivery_cost_per_FVP = ifelse(grepl("-", vacc_delivery_cost_per_FVP), sum(as.numeric(unlist(strsplit(vacc_delivery_cost_per_FVP, split="-"))))/2, as.numeric(vacc_delivery_cost_per_FVP))) %>%
      dplyr::mutate(orig_delivery_cost_per_dose = ifelse(grepl("-", vacc_delivery_cost_per_dose), sum(as.numeric(unlist(strsplit(vacc_delivery_cost_per_dose, split="-"))))/2, as.numeric(vacc_delivery_cost_per_dose))) %>%
      dplyr::mutate(orig_procurement_cost_per_FVP = ifelse(grepl("-", vacc_procurement_cost_per_FVP), sum(as.numeric(unlist(strsplit(vacc_procurement_cost_per_FVP, split="-"))))/2, as.numeric(vacc_procurement_cost_per_FVP))) %>%
      dplyr::mutate(orig_procurement_cost_per_dose = ifelse(grepl("-", vacc_procurement_cost_per_dose), sum(as.numeric(unlist(strsplit(vacc_procurement_cost_per_dose, split="-"))))/2, as.numeric(vacc_procurement_cost_per_dose))) %>%
      dplyr::mutate(vacc_delivery_cost_per_FVP = get_inflation_multiplier(USD_year)*orig_delivery_cost_per_FVP) %>%
      dplyr::mutate(vacc_delivery_cost_per_dose = get_inflation_multiplier(USD_year)*orig_delivery_cost_per_dose) %>%
      dplyr::mutate(vacc_procurement_cost_per_FVP = get_inflation_multiplier(USD_year)*orig_procurement_cost_per_FVP) %>%
      dplyr::mutate(vacc_procurement_cost_per_dose = get_inflation_multiplier(USD_year)*orig_procurement_cost_per_dose) %>%
      dplyr::mutate(vacc_total_cost_per_dose = vacc_delivery_cost_per_dose + vacc_procurement_cost_per_dose) %>%
      dplyr::mutate(vacc_total_cost_per_FVP = vacc_delivery_cost_per_FVP + vacc_procurement_cost_per_FVP) %>%
      ungroup

    if (sensitivity == "base"){
      delivery_cost_per_FVP <- median(cost_dat$vacc_delivery_cost_per_FVP, na.rm = TRUE)
      procurement_cost_per_FVP <- median(cost_dat$vacc_procurement_cost_per_FVP, na.rm = TRUE)
      total_cost_per_FVP <- median(cost_dat$vacc_total_cost_per_FVP, na.rm = TRUE)
    } else if (sensitivity == "low"){
      delivery_cost_per_FVP <- unname(stats::quantile(cost_dat$vacc_delivery_cost_per_FVP, probs = c(0.025), na.rm=TRUE))
      procurement_cost_per_FVP <- unname(stats::quantile(cost_dat$vacc_procurement_cost_per_FVP, probs = c(0.025), na.rm=TRUE))
      total_cost_per_FVP <- unname(stats::quantile(cost_dat$vacc_total_cost_per_FVP, probs = c(0.025), na.rm=TRUE))
    } else if (sensitivity == "high"){
      delivery_cost_per_FVP <- unname(stats::quantile(cost_dat$vacc_delivery_cost_per_FVP, probs = c(0.975), na.rm=TRUE))
      procurement_cost_per_FVP <- unname(stats::quantile(cost_dat$vacc_procurement_cost_per_FVP, probs = c(0.975), na.rm=TRUE))
      total_cost_per_FVP <- unname(stats::quantile(cost_dat$vacc_total_cost_per_FVP, probs = c(0.975), na.rm=TRUE))
    }   

    cost_ls <- list(sensitivity = sensitivity, deliv_per_FVP = delivery_cost_per_FVP, proc_per_FVP = procurement_cost_per_FVP, tot_per_FVP = total_cost_per_FVP)
  return(cost_ls)
}


#' @name get_vacc_cost_df
#' @title get_vacc_cost_df
#' @description Create dataframe to get vaccination costs per FVP for baseline, lower and upper bound sensitivity parameters
#' @param data_wd
#' @param sensitivity
#' @return 
get_vacc_cost_df <- function(){
  fvp_cost_base <- get_vacc_cost_data(sensitivity = "base")
  fvp_cost_low <- get_vacc_cost_data(sensitivity = "low")
  fvp_cost_high <- get_vacc_cost_data(sensitivity = "high")
  bind_df <- bind_rows(unlist(fvp_cost_base), unlist(fvp_cost_low), unlist(fvp_cost_high))
  fvp_cost_df <- bind_df %>% 
    tidyr::gather(cost_type, cost_per_fvp, contains("per_FVP")) %>%
    rowwise %>% 
    dplyr::mutate(cost_type = gsub("_per_FVP", "", cost_type)) %>% ungroup %>%
    dplyr::mutate(cost_per_fvp = as.numeric(cost_per_fvp))

  return(fvp_cost_df)
}

#' @name get_inflation_multiplier
#' @title get_inflation_multiplier
#' @description Adjust for inflation to 2017 USD dollars according to World Bank CPI
#' @param data_wd
#' @param sensitivity
#' @return
get_inflation_multiplier <- function(in_year){
  
  us_cpi <- data.frame(year = 2007:2016, cpi = c(95.08699238, 98.73747739, 98.38641997, 100, 103.1568416, 105.2915045, 106.8338489, 108.5669321, 108.695722, 110.0670089)) ## data from World Bank Consumer Price Index
  cpi2017 <- 112.4115573
  us_cpi <- us_cpi %>%
    dplyr::mutate(cum_inf_rate_multiplier = 1+((cpi2017-cpi)/cpi))

  multiplier <- us_cpi %>% 
    dplyr::filter(year == in_year) %>% 
    dplyr::select(cum_inf_rate_multiplier) %>% unname %>% unlist

  return(multiplier) 
}


#' @name get_discounted_cost_df
#' @title get_discounted_cost_df
#' @description Get monetary costs discounted by 3% in 2017 dollars.
#' @param current_cost
#' @return
get_discounted_cost_df <- function(current_cost){
  
  cost_adj_df <- map_dfr(2018:2035, function(year){
    cost_adj = current_cost/((1+0.03)^(year-2018))
    data.frame(year = year, cost = cost_adj)
  })

  return(cost_adj_df) 
}


#' @name get_country_gdp
#' @title get_country_gdp
#' @description Import WorldBank 2017 GDP data in 2017 dollars.
#' @param infile
#' @return
get_country_gdp <- function(infile = "data/WorldBank_GDP/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_10473720.csv"){
  
  wbdat <- read_csv(infile, skip = 5, col_types = paste0("_c__", paste0(rep("_",57), collapse=""),"d__"), col_names = c("cntry_code", "gdp2017"))

  return(wbdat) 
}


#' @name get_vacc_coverage
#' @title get_vacc_coverage
#' @description Helper function to read vaccination campaign coverage estimates and average relevant values (2-dose estimates). Estimates are either unweighted or weighted by the total number of vaccinated individuals contributing to a coverage estimate.
#' @param data_file Spreadsheet where coverage survey data were reviewed
#' @return dataframe with estimates for unweighted and weighted mean vaccination campaign coverage proportion across coverage surveys in the review
get_vacc_coverage <- function(data_file = "data/Review_OCV_campaign_coverage.ods"){

  set.seed(35546660)

  covDat <- read_ods(data_file, sheet = "coverage_survey", col_names = TRUE, na = "NA") %>%
    dplyr::filter(redundant_aggregation == FALSE & two_doses == TRUE & !exclude)
  ci_low <- sapply(strsplit(covDat$vacc_coverage_ci, split = "-"), function(x){return(x[1])})
  ci_high <- sapply(strsplit(covDat$vacc_coverage_ci, split = "-"), function(x){return(x[2])})

  clDat <- covDat  %>%
    dplyr::mutate(sample_size = round(number_vacc/(vacc_coverage_pct/100))) %>%
    dplyr::mutate(vacc_coverage = vacc_coverage_pct/100, vacc_coverage_ci_lower = as.numeric(ci_low)/100, vacc_coverage_ci_upper = as.numeric(ci_high)/100) %>%
    dplyr::mutate(vacc_coverage_inverseVar = 1/((sqrt(sample_size)*(vacc_coverage_ci_upper-vacc_coverage_ci_lower)/3.92)^2)) %>%
    dplyr::mutate(coverage_se_norm = (vacc_coverage_ci_upper-vacc_coverage_ci_lower)/(1.96*2))
  
  distDat <- clDat %>% 
    dplyr::select(study_id, vaccinated_group, vacc_coverage, coverage_se_norm) 

  samps <- unlist(unlist(lapply(unique(distDat$study_id), function(id){
    studyDat <- distDat %>% dplyr::filter(study_id==id)
    sampleStudyDat <- sample_n(studyDat, size = 5000, replace = TRUE)
    vals <- unlist(lapply(1:nrow(sampleStudyDat), function(i){
      row <- unlist(sampleStudyDat[i,3:4])
      rnorm(n = 1, mean = row['vacc_coverage'], sd = row['coverage_se_norm'])
    })) 
    return(vals)
  })))

  ## wMedian1090 was used for sensitivity analysis
  med1090_ci_vec <- c("version" = "wMedian1090", 
                  "middle" = median(samps), 
                  "ci_lower" = unname(quantile(samps, probs=.1)), 
                  "ci_upper" = unname(quantile(samps, probs=.9)))
  mnIQR_ci_vec <- c("version" = "wMeanIQR", 
                  "middle" = mean(samps), 
                  "ci_lower" = unname(quantile(samps, probs=.25)), 
                  "ci_upper" = unname(quantile(samps, probs=.75)))
  medIQR_ci_vec <- c("version" = "wMedianIQR", 
                  "middle" = median(samps), 
                  "ci_lower" = unname(quantile(samps, probs=.25)), 
                  "ci_upper" = unname(quantile(samps, probs=.75)))

  summDat <- rbind(mnIQR_ci_vec, medIQR_ci_vec, med1090_ci_vec)
  summDat <- tbl_df(summDat) %>%
    tidyr::gather(quant, estimate, middle:ci_upper) %>%
    dplyr::mutate(method = paste(version, quant, sep = " "))

  return(summDat)
}


#' @name generate_pct_protect_function
#' @title generate_pct_protect_function
#' @description Using vaccine efficacy studies in Bi et al. (2017), generate a function that takes the year since vaccination and provides an estimate of the direct ve. Vaccine efficacy declines to 0 after my_trunc_year years
#' @param my_trun_year  
#' @param my_ve_scen 
#' @return dataframe with estimates for unweighted and weighted mean vaccination campaign coverage proportion across coverage surveys in the review
generate_pct_protect_function <- function(my_trunc_year = 5, my_ve_scen = "base"){

    ve.dat <- read_csv("data/ocv_ve_overtime.csv")

    ve.dat$T <- (ve.dat$TL+ve.dat$TR)/2

    ve.dat$weights <-  1/(abs(ve.dat$se)^2)

    ve.dat <- ve.dat[ve.dat$TL<48,]

    ##this is our basic trend in vaccine effictiveness.
    ve.trend <- lm(yi~T , data=ve.dat, weights = weights)

    ##create a function that gives the expected percent protected
    ##by a vaccine during a particular year after vaccination

    pct.protect<-function(year, ci=FALSE,trunc_year=my_trunc_year,ve_scen=my_ve_scen) {
        if (any(year%%1 !=0)) {warning("function designed to average across years only")}
        if(!ve_scen %in% c("base","low","high")) {warning("ve_scen not recognized, assuming ve_scen='base'")}

        months <- (year-.5)*12

        if(ve_scen %in% c("high","low")){
            ci <- TRUE
        }

        if (ci) {
            rc <- pmax(1-exp(predict(ve.trend, newdata=data.frame(T=months), interval="confidence") ),0) %>% unname
        } else {
            rc <- pmax(1-exp(predict(ve.trend, newdata=data.frame(T=months))),0) %>% unname
        }

        rc <- as.matrix(rc)

        if(any(year>trunc_year)){
            rc[which(year>trunc_year),] <- 0
        }

        if(ve_scen == "high"){
            rc <- rc[,2]
        } else if (ve_scen == "low"){
            rc <- rc[,3]
        }

        return(rc)
    }
    return(pct.protect)
}



#' @name generate_indirect_incidence_mult
#' @title generate_indirect_incidence_mult
#' @description Using studies on indirect vaccine protection from Kolkota and Matlab, generate a function that takes the level of vaccination coverage and returns a multiplier indicating the percentage reduction in incidence due to indirect vaccine protection. These represent the baseline parameters for indirect vaccine protection. 
#' @return 
generate_indirect_incidence_mult <- function(){

    ## data from Kolkota trial
    vc_u_k <- c(25,28,31,34,45) # upper limit of coverage bins
    vc_l_k <- c(0,25,28,31,34) # lower limit of coverage bins
    risk_v_k <- c(1.28,1.48,1.01,0.97,1.24) # incid among vacc recipients
    risk_p_k <- c(5.54,5.64,2.48,2.25,1.93) # incid among placebo recipients
    
    ## data from Matlab trial
    vc_u_m <- c(28,35,40,50,60) # upper limit of coverage bins
    vc_l_m <- c(0,28,35,40,50) # lower limit of coverage bins
    risk_v_m <- c(2.66,2.47,1.57,2.25,1.27) # incid among vacc recipients
    risk_p_m <- c(7.01,5.87,4.72,4.65,1.47) # incid among placebo recipients

                                        #effective_coverage <- function(cov){ cov*.75 + .25}

    df = data.frame(coverage=c((vc_u_k+vc_l_k)/2, (vc_u_m+vc_l_m)/2)/100,
                    vaccinated=c(risk_v_k, risk_v_m),
                    placebo=c(risk_p_k, risk_p_m),
                    ve_I=c(1-risk_v_k/risk_p_k, 1-risk_v_m/risk_p_m),
                    loc=c(rep("Kolkata", length(vc_u_k)), rep("Matlab", length(vc_u_m))))

    df <- df %>% group_by(loc) %>%
      dplyr::mutate(indirect=pmax(0.001, 1 - (placebo/placebo[1])), effective_cov=coverage) %>% 
      ungroup

    ##df %>%  ggplot(aes(x=coverage,y=indirect)) + geom_point(aes(color=loc))

    fit <- lm(log(indirect/(1-indirect)) ~ effective_cov, data=df)
    pred <- predict(fit, newdata = data.frame(effective_cov=seq(0.01,1,length=100)))

    indirect_inc_mult <- function(effective_coverage){
        tmp <- predict(fit, newdata = data.frame(effective_cov=effective_coverage))
        rc <- (1- (exp(tmp)/(1+exp(tmp)))) %>% unname

        rc[which(effective_coverage==0)] <- 1
        rc[which(effective_coverage==1)] <- 0
        # returns a multiplier for incidence based on coverage 
        # (if effective coverage is 0%, there is no reduction due to indirect effects and multiplier is 1)
        return(rc)
    }

}

#' @name generate_indirect_incidence_sensitivity_mult
#' @title generate_indirect_incidence_sensitivity_mult
#' @description Generate a function that takes the level of vaccination coverage and returns a multiplier indicating the percentage reduction in incidence due to indirect vaccine protection. This function is for the two sensitivity analyses. Lower bound analyses have no indirect protection and upper bound analyses follow Longini Jr., et al. (2007) data.
#' @param bound upper or lower 
#' @return 
generate_indirect_incidence_sensitivity_mult <- function(bound){

  if(bound == "lower"){
    indirect_inc_mult <- function(effective_coverage){
        ## lower bound is no indirect effects for any cell
        return(1)
    }
  } else if(bound == "upper"){
    ## Use indirect effects estimates from "Controlling Endemic Cholera with Oral Vaccines", Longini Jr., et al. (2007) PLOS MED 
    ## See Table 2 mean indirect effectiveness
    df <- data.frame(effective_cov=c(10,30,50,70,90)/100, indirect=c(.3,.7,.89,.97,.99))
    fit <- lm(log(indirect/(1-indirect)) ~ effective_cov, data=df)
    pred <- predict(fit, newdata = data.frame(effective_cov=seq(0.01,1,length=100)))

    indirect_inc_mult <- function(effective_coverage){
      tmp <- predict(fit, newdata = data.frame(effective_cov=effective_coverage))
      rc <- (1- (exp(tmp)/(1+exp(tmp)))) %>% unname

      rc[which(effective_coverage==0)] <- 1
      rc[which(effective_coverage==1)] <- 0
      # returns a multiplier for incidence based on coverage 
      # (if effective coverage is 0%, there is no reduction due to indirect effects and multiplier is 1)
      return(rc)
    }
  }

  return(indirect_inc_mult)
}


#' @name generate_flatline_multiplier
#' @title generate_flatline_multiplier
#' @description Generate a function that represents projected secular trends in cholera incidence. This multiplier is used to adjust projected cholera incidence in future years.
#' @return
generate_flatline_multiplier <- function(){
  flatline_multiplier <- function(year, base_year = 2016) { return(1) }
  return(flatline_multiplier)
}


#' @name load_life_expect
#' @title load_life_expect
#' @description Import and process 2017 UN World Population Prospects estimates on projected life expectancy
#' @param file
#' @return
load_life_expect <- function(file="data/WPP2017_SA4_MORT_F07_1_LIFE_EXPECTANCY_0_BOTH_SEXES.xlsx"
                             ){

    dat <- read_xlsx(file,skip = 16) %>%
    rename(country="Type of aggregate, group, and constituents *") %>%
    dplyr::select(country,`1990-1995`,`1995-2000`,`2000-2005`,`2005-2010`,`2010-2015`,`2015-2020`,`2020-2025`,`2025-2030`,`2030-2035`,`2035-2040`)

    le_dat <- bind_rows(
        dat %>%
        tidyr::gather(year,life_expect,-country) %>%
        dplyr::mutate(year=substr(year,1,4) %>% as.numeric),
        dat %>%
        tidyr::gather(year,life_expect,-country) %>%
        dplyr::mutate(year=(substr(year,1,4) %>% as.numeric)+1),
        dat %>%
        tidyr::gather(year,life_expect,-country) %>%
        dplyr::mutate(year=(substr(year,1,4) %>% as.numeric)+2),
        dat %>%
        tidyr::gather(year,life_expect,-country) %>%
        dplyr::mutate(year=(substr(year,1,4) %>% as.numeric)+3),
        dat %>%
        tidyr::gather(year,life_expect,-country) %>%
        dplyr::mutate(year=(substr(year,1,4) %>% as.numeric)+4))

    ## lookup table doesn't seem to be working perfectly
    le_dat <- le_dat %>% mutate(country=recode(country,"C\u00f4te d'Ivoire"="Cote d'Ivoire"))

    ## match country names and codes and drop blanks
    country_names <- unique(le_dat$country)
    fixed_names <- data.frame(country=country_names,country_code=sapply(country_names,fix_country_name,verbose=TRUE) %>% unlist) %>%
    dplyr::mutate(
               country_code=ifelse(nchar(as.character(country_code))==3,as.character(country_code),NA))

    le_dat <- left_join(le_dat,fixed_names) %>% filter(!is.na(country_code))

    return(le_dat)
}


#' @name load_life_expect
#' @title load_life_expect
#' @description Return life expectancy for a specific country and year
#' @param cntry vector of countries
#' @param year  vector of years
#' @return
get_life_expect <- function(cntry, year){
    if(length(cntry) ==1 & length(year)>1){
        cntry <- rep(cntry,length(year))
    } else if(length(year) == 1 & length(cntry)>1){
        year <- rep(year,length(cntry))
    } else if(length(year) != length(cntry)){
        stop("length of year and cntry should match, or one should be length 1")
    }

    cat(sprintf("loading life expectancy table \n"))
    le_tab <- load_life_expect() %>% distinct()

    #recover()
    rc <- left_join(data.frame(country_code=sapply(cntry,fix_country_name) %>% as.character,
                               year=round(year)),
                    le_tab) %>%
    mutate(birth_year=year) %>%
    dplyr::select(country_code,birth_year,life_expect) %>%
      rename(cntry_code=country_code)
    return(rc)
}


##' gets average age of infection for a given country
##' @param cntry
##' @return
get_inf_age <- function(cntry){
    cntry_code <- sapply(cntry,fix_country_name)

    ## values for average age derived from this code:
    ## library(metafor)
    ## age_cases <- read_xlsx("data/age_cases.xlsx")

    ## ## using relationship between variacne and mean of age to get variance of log(age)
    ## ## 10.1002/sim.1525
    ## age_cases %>% filter(!is.na(age_sd)) %>%
    ## mutate(mean_log=log(mean_age),var_log=log(1+age_sd^2/mean_age^2)) %>%
    ## summarize(mean_age_inf=weighted.mean(mean_age,w=1/var_log)) %>% unlist

    ## aoi_south_asia <-  age_cases %>% filter(!is.na(age_sd),location %in% c('IND','BGD')) %>%
    ## mutate(mean_log=log(mean_age),var_log=log(1+age_sd^2/mean_age^2)) %>%
    ## summarize(mean_age_inf=weighted.mean(mean_age,w=1/var_log)) %>% unlist

    ## aoi_africa_haiti <- age_cases %>% filter(!is.na(age_sd),!location %in% c('IND','BGD')) %>%
    ## mutate(mean_log=log(mean_age),var_log=log(1+age_sd^2/mean_age^2)) %>%
    ## summarize(mean_age_inf=weighted.mean(mean_age,w=1/var_log)) %>% unlist


    inf_age <- case_when(
        cntry_code %in% c("BGD","IND","NPL") ~ 22.42,
        !cntry_code %in% c("BGD","IND","NPL") ~ 25.75
    )

    return(inf_age)
}


##' gets assumed cfr for a given country for all countries in study area
##' @param cntry
##' @regen_who_ests
##' @return
get_cfr <- function(cntry){

  cntry_code <- sapply(cntry, fix_country_name) %>% data.frame
  colnames(cntry_code) <- "ISO_A1"
  cntry_code[,1] <- as.character(cntry_code[,1])

  deaths_summary <- read_csv("data/who_cfrs.csv")

  total_cfrs <- deaths_summary %>% 
    dplyr::filter(!is.na(cases), !is.na(deaths), cases>0) %>% 
    dplyr::rename(cfr_tot = cfr) %>%
    dplyr::select(cntry_code, cfr_tot)

  total_cfrs_select <- total_cfrs %>%
    dplyr::rename(ISO_A1 = cntry_code) %>%
    dplyr::mutate(cfr_tot = ifelse(cfr_tot>=.07, NA, cfr_tot)) 
  mn_cfr <- mean(total_cfrs_select$cfr_tot, na.rm = TRUE)

  rc <- left_join(cntry_code, total_cfrs_select, by = c("ISO_A1")) %>%
    dplyr::mutate(cfr = ifelse(is.na(cfr_tot), mn_cfr, cfr_tot)) %>%
    dplyr::select(ISO_A1, cfr)

  return(rc)

}



##' summarizes number of vaccines allocated to districts
##' @param targets_dir
##' @return
get_target_district_summary <- function(targets_dir = "country_targets"){

  ## 1/1/2019 rm cell strategy code
  fns <- list.files(targets_dir)
  districtTargets <- paste0(targets_dir,"/",fns[grep("districtTargets.csv$",fns,perl=T)])
  cellTargets <- paste0(targets_dir,"/",fns[grep("cellTargets.csv$",fns,perl=T)])
  ## cell targets for district targeting allocation strategies
  cellTargets_for_districtStrategy <- cellTargets[grepl("District", cellTargets)]

  if(length(cellTargets_for_districtStrategy)>0){
    cellSummD <- map_dfr(cellTargets_for_districtStrategy, function(x) read_csv(x) %>% 
                                          dplyr::mutate(scenario=gsub("_[0-9]{4}_cellTargets.csv", "", gsub(paste0(targets_dir, "/"), "", x)))) %>%
              dplyr::mutate(uqid = paste(scenario, district, year, sep = "_")) %>%
              dplyr::filter(fvp > 0) %>%
              group_by(uqid) %>%
              summarise(vacc_used = sum(vacc_used, na.rm = TRUE), 
                        fvp = sum(fvp, na.rm = TRUE)) %>%
              ungroup 
  } 

  if(length(districtTargets)>0){
    districtSumm <- map_dfr(districtTargets,function(x) read_csv(x) %>% 
                                                    dplyr::mutate(scenario=substr(gsub("_districtTargets.csv", "", gsub(paste0(targets_dir, "/"), "", x)), 1, nchar(gsub("_districtTargets.csv", "", gsub(paste0(targets_dir, "/"), "", x)))-5))) %>%
                    dplyr::mutate(uqid = paste(scenario, district, year, sep = "_")) %>%
                    dplyr::select(scenario, cntry_code, year, district, incidence, cases, pop, index, id, vacc_used, fvp, uqid)
  }

  ## add values from summarized cell targets to appropriate district-year
  for (i in 1:nrow(cellSummD)){
    uqid <- cellSummD[i,]$uqid
    vacc_used <- cellSummD[i,]$vacc_used
    fvp <- cellSummD[i,]$fvp
    districtSumm[which(districtSumm$uqid == uqid),]$vacc_used <- vacc_used
    districtSumm[which(districtSumm$uqid == uqid),]$fvp <- fvp
  }

  targets_df <- districtSumm %>% dplyr::filter(fvp > 0)
  return(targets_df)
}


##' summarizes number of vaccines allocated to countries
##' @param targets_dir
##' @return
get_target_summary <- function(targets_dir="country_targets"){

  fns <- list.files(targets_dir)
  districtTargets <- paste0(targets_dir,"/",fns[grep("districtTargets.csv$",fns,perl=T)])
  cellTargets <- paste0(targets_dir,"/",fns[grep("cellTargets.csv$",fns,perl=T)])
  ## cell targets for district targeting allocation strategies
  cellTargets_for_districtStrategy <- cellTargets[grepl("District", cellTargets)]
  ## cell targets for cell targeting allocation strategies
  cellTargets_for_cellStrategy <- cellTargets[grepl("Cell", cellTargets)]

  if(length(cellTargets_for_districtStrategy)>0){
    cellSummD <- map_dfr(cellTargets_for_districtStrategy, function(x) read_csv(x) %>% 
                                          dplyr::mutate(scenario=gsub("_[0-9]{4}_cellTargets.csv", "", gsub(paste0(targets_dir, "/"), "", x)))) %>%
              dplyr::mutate(uqid = paste(scenario, cntry_code, year, sep = "_")) %>%
              group_by(uqid) %>%
              summarise(scenario = first(scenario),
                        cntry_code = first(cntry_code),
                        year = first(year), 
                        vacc_used = sum(vacc_used, na.rm = TRUE),
                        fvp = sum(fvp, na.rm = TRUE), 
                        target = "partial") %>%
              ungroup 
  } 
  

  if(length(cellTargets_for_cellStrategy)>0){
    cellSummC <- map_dfr(cellTargets_for_cellStrategy, function(x) read_csv(x) %>%
                                            dplyr::mutate(scenario=gsub("_[0-9]{4}_cellTargets.csv", "", gsub(paste0(targets_dir, "/"), "", x)))) %>%
              dplyr::mutate(uqid = paste(scenario, cntry_code, year, sep = "_")) %>%
              group_by(uqid) %>%
              summarise(scenario = first(scenario),
                        cntry_code = first(cntry_code),
                        year = first(year),
                        vacc_used = sum(vacc_used, na.rm = TRUE),
                        fvp = sum(fvp, na.rm = TRUE)) %>%
              ungroup
  } 

  if(length(districtTargets)>0){
    districtSumm <- map_dfr(districtTargets,function(x) read_csv(x) %>% 
                                                    dplyr::mutate(scenario=substr(gsub("_districtTargets.csv", "", gsub(paste0(targets_dir, "/"), "", x)), 1, nchar(gsub("_districtTargets.csv", "", gsub(paste0(targets_dir, "/"), "", x)))-5))) %>%
                    dplyr::mutate(uqid = paste(scenario, district, year, sep = "_")) %>%
                    dplyr::mutate(target = ifelse(fullDistrict_target == TRUE, "complete", NA)) %>%
                    dplyr::select(scenario, cntry_code, year, district, incidence, cases, pop, index, id, vacc_used, fvp, target, uqid)
  } 
   

    ####################################################
    #### number of targets and fvp per scenario ####
    ## district scenarios ##
    numTargetsDd <- districtSumm %>% 
      dplyr::filter(vacc_used > 0) %>%
      group_by(scenario, cntry_code, year) %>%
      summarise(numTargets = n(), ddFVP = sum(fvp)) %>%
      ungroup
    numTargetsDc <- cellSummD %>%
      dplyr::filter(vacc_used > 0) %>%
      group_by(scenario, cntry_code, year) %>%
      summarise(dcFVP = sum(fvp)) %>%
      ungroup
    numTargetsD <- full_join(numTargetsDd, numTargetsDc, by = c("scenario", "cntry_code", "year")) %>%
      dplyr::mutate(dcFVP = ifelse(is.na(dcFVP), 0, dcFVP), ddFVP = ifelse(is.na(ddFVP), 0, ddFVP)) %>%
      dplyr::mutate(fvp = dcFVP + ddFVP) %>%
      dplyr::select(scenario, cntry_code, year, numTargets, fvp)

    ## cell scenarios ##
    if(exists("cellSummC")){
      numTargetsC <- cellSummC %>%
      dplyr::filter(vacc_used > 0) %>%
      group_by(scenario, cntry_code, year) %>%
      summarise(numTargets = nrow(.), fvp = sum(fvp)) 
      numTargets_df <- bind_rows(numTargetsD, numTargetsC)
    } else{
      numTargets_df <- numTargetsD
    }
    

    return(numTargets_df)

}


##' summarizes number of vaccines allocated to each district for a district-level strategy
##' @param targets_dir
##' @return
get_district_targets_districtStrategy <- function(targets_dir="country_targets"){

  fns <- list.files(targets_dir)
  districtTargets <- paste0(targets_dir,"/",fns[grep("districtTargets.csv$",fns,perl=T)])
  cellTargets <- paste0(targets_dir,"/",fns[grep("cellTargets.csv$",fns,perl=T)])
  
  ## aggregate cell targets for district targeting allocation strategies to district level
  cellTargets_for_districtStrategy <- cellTargets[grepl("District", cellTargets)]
  cellSummD <- map_dfr(cellTargets_for_districtStrategy, function(x) read_csv(x) %>% 
                                          dplyr::mutate(scenario=gsub("_[0-9]{4}_cellTargets.csv", "", gsub(paste0(targets_dir, "/"), "", x)))) %>%
              group_by(scenario, district, year) %>%
              summarise(vacc_used = sum(vacc_used, na.rm = TRUE), fvp = sum(fvp, na.rm = TRUE), target = "partial") %>%
              dplyr::filter(vacc_used > 0) %>%
              dplyr::mutate(uqid = paste(scenario, district, year, sep = "_")) %>%
              dplyr::select(uqid, vacc_used, fvp, target) %>%
              ungroup

  districtSumm <- map_dfr(districtTargets,function(x) read_csv(x) %>% 
                                                    dplyr::mutate(scenario=substr(gsub("_districtTargets.csv", "", gsub(paste0(targets_dir, "/"), "", x)), 1, nchar(gsub("_districtTargets.csv", "", gsub(paste0(targets_dir, "/"), "", x)))-5))) %>%
                  dplyr::mutate(uqid = paste(scenario, district, year, sep = "_")) %>%
                  dplyr::mutate(target = ifelse(fullDistrict_target == TRUE, "complete", NA)) %>%
                  dplyr::select(scenario, cntry_code, year, district, incidence, cases, pop, index, id, vacc_used, fvp, target, uqid)

  for (i in 1:nrow(cellSummD)){
    districtSumm[which(districtSumm$uqid == cellSummD$uqid[i]),c("vacc_used", "fvp", "target")] <- cellSummD[i, c("vacc_used", "fvp", "target")]
  }

  return(districtSumm)
  
}


##' generates gathered data frame of raw samples by country, year and scenario
##' @param scenarioLs
##' @param gen_output_dir
##' @param cf_wd
##' @return
get_run_samples <- function(scenarioLs, 
                            gen_output_dir="generated_outputs/",
                            cf_wd = "generated_outputs/cf/"){
  
  dummy <- list.files(path = cf_wd, pattern = "_cf_cases_samples")
  countryLs <- sort(unique(substring(dummy, 1, 3)))

  ##loop countries
  country_df <- map_dfr(countryLs, function(cntry){

    print(paste("processing", cntry, "samples. please be patient..."))
    ## import health outcome data
    cfr <- get_cfr(cntry)[,2]
    aoi <- get_inf_age(cntry)
    dis_weight=0.247
    inf_dur=4/365

    ## loop scenarios
    scens_df <- map_dfr(scenarioLs, function(scenario){

      print(paste("    ", scenario, "scenario"))
      alloc_strategy <- unlist(strsplit(scenario, split="_"))[1]

      case_fn <- list.files(path = gen_output_dir, pattern = paste0(cntry, "_", scenario, "_cases_samples"), full.names = TRUE)
      cf_fn <- list.files(path = cf_wd, pattern = paste0(cntry, "_cf_cases_samples"), full.names = TRUE)
		print(case_fn)
		print(cf_fn)
      ## calculate health outcomes for all years
      case_samples <- read_csv(case_fn) %>%
        dplyr::select(starts_with("x20"))
		print(head(case_samples))
      cf_samples <- read_csv(cf_fn) %>%
        dplyr::select(starts_with("x20"))
		print(head(cf_samples))
      casesAverted_samples <- cf_samples - case_samples
      var_years <- names(casesAverted_samples)

      ## loop years
      years_df <- map_dfr(1:length(var_years), function(i){

        print(paste("        ", var_years[i]))

        yr <- as.numeric(substring(var_years[i], 2, 5))
        ## year-specific health outcome data
        birth_year <- round(yr-aoi)
        life_expect_df <- get_life_expect(cntry, birth_year) %>% distinct
        ## set average age of infection to lower value between life expectancy and avg age of infection table (RWA genocide)
        aoi_cl <- pmin(get_inf_age(cntry), get_life_expect(cntry, birth_year)$life_expect) 

        ## calculate outcomes for all posterior samples
        case_yr <- case_samples %>% dplyr::select(!!var_years[i]) %>% unlist
        cf_yr <- cf_samples %>% dplyr::select(!!var_years[i]) %>% unlist
        avert_yr <- casesAverted_samples %>% dplyr::select(!!var_years[i]) %>% unlist
        deaths_averted_summ <- avert_yr * cfr
        yll_averted_summ <- deaths_averted_summ * (life_expect_df$life_expect - aoi_cl)
        yld_averted_summ <- avert_yr * inf_dur * dis_weight
        dalys_averted_summ <- yll_averted_summ + yld_averted_summ

        gather_df <- data.frame(cntry_code = cntry, 
                          scenario = scenario,
                          alloc_strategy = alloc_strategy,
                          year = yr,
                          samp = 1:length(avert_yr),
                          cases = case_yr,
                          cases_cf = cf_yr, 
                          casesAverted = avert_yr,
                          deathsAverted = deaths_averted_summ,
                          dalysAverted = dalys_averted_summ) %>%
          tidyr::gather(measure, value, cases:dalysAverted)
        return(gather_df)
      })
      return(years_df)
    })
    return(scens_df)
  })
  return(country_df)
}


##' generates gathered data frame of cumulative samples by country, year and scenario
##' @param import_df output from get_run_samples
##' @return
process_cum_run_samples <- function(import_df){

  print("processing cumulative scenario samples. please be patient.")
  prep_cum_df <- spread(import_df %>% dplyr::mutate(year = paste0("x", year)), year, value)
  rm(import_df)
  gc()

  metacol <- prep_cum_df %>% dplyr::select(cntry_code, scenario, alloc_strategy, measure, samp)
  dataonly <- prep_cum_df %>% dplyr::select(contains("x"))
  cum_df <- data.frame(bind_cols(data.frame(t(apply(dataonly, 1, cumsum))), metacol)) %>%
    tidyr::gather(year, value, contains("x")) %>%
    dplyr::mutate(year = as.numeric(gsub("x", "", year)))

  rm(dataonly)
  gc()

  return(cum_df)
}

##' prepare median and CI for entire study area from processed raw or cumulative samples
##' @param processed_samples_df
##' @param scenarioLs
##' @return
process_CI_ssa <- function(processed_samples_df, scenarioLs){

  ssa_samples <- processed_samples_df %>%
    dplyr::mutate(scenario = as.character(scenario)) %>%
    dplyr::filter(scenario %in% scenarioLs) %>%
    dplyr::select(-cntry_code) %>%
    group_by(scenario, measure, year, samp) %>%
    dplyr::summarise(value = sum(value)) 
  rm(processed_samples_df)
  gc()

  ## process raw samples
  ssa_samples_med <- ssa_samples %>% 
    group_by(scenario, measure, year) %>% 
    dplyr::summarise(median = median(value), mean = mean(value)) %>% 
    ungroup
  ssa_samples_lower <- ssa_samples %>% 
    group_by(scenario, measure, year) %>% 
    dplyr::summarise(ci_lower = quantile(value, probs = 0.025)) %>% 
    ungroup
  ssa_samples_upper <- ssa_samples %>% 
    group_by(scenario, measure, year) %>% 
    dplyr::summarise(ci_upper = quantile(value, probs = 0.975))%>% 
    ungroup
  rm(ssa_samples)
  gc()

  ssa_samples_summ <- full_join(ssa_samples_med, ssa_samples_lower, by = c("scenario", "measure", "year")) %>%
    full_join(ssa_samples_upper, by = c("scenario", "measure", "year")) %>%
    rowwise %>%
    dplyr::mutate(alloc_strategy = unlist(strsplit(scenario, split = "_"))[1]) %>%
    ungroup

  return(ssa_samples_summ)

}


##' generate wash strategy labels based on scenario
##' @return
get_washname <- function(scenario) {
  ifelse(grepl("San", scenario), "san", ifelse(grepl("Wat", scenario), "wat", ifelse(grepl("Wash", scenario), "wash", NA)))
}


##' generate health outcome labels for plotting
##' @return
label_healthOutcomes <- function(){

  label_df <- data.frame(measure = c("cases", "cases_cf", "casesAverted", "deathsAverted", "yllAverted", "yldAverted", "dalysAverted", "fvp", "caPer1kFVP", "hr_pop"), 
                        plotMeasure = c("Cases", "Cases (no vacc.)", "Cases Averted", "Deaths Averted", "YLL Averted", "YLD Averted", "DALYs Averted", "Fully vacc. persons", "Cases Averted per 1000 FVP", "High-Risk Population"))

  return(label_df)
}


##' generate country name labels for plotting
##' @return
label_countries <- function(){

  country_labs <- data.frame(cntry_code = c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GAB", "GHA", "GIN", "GMB", "GNB", "GNQ", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SDN", "SEN", "SLE", "SOM", "SSD", "SWZ", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE"), country_name = c("Angola", "Burundi", "Benin", "Burkina Faso", "C Afr Rep", "Cote D'Ivoire", "Cameroon", "DR Congo", "Rep of Congo", "Ethiopia", "Gabon", "Ghana", "Guinea", "Gambia", "Guinea-Bissau", "Eq Guinea", "Kenya", "Liberia", "Madagascar", "Mali", "Mozambique", "Mauritania", "Malawi", "Namibia", "Niger", "Nigeria", "Rwanda", "Sudan", "Senegal", "Sierra Leone", "Somalia", "S Sudan", "Swaziland", "Chad", "Togo", "Tanzania", "Uganda", "S Africa", "Zambia", "Zimbabwe"), stringsAsFactors = FALSE)

  return(country_labs)
}


##' generate allocation strategy labels for plotting
##' @return
label_allocStrategy <- function(){

  label_df <- data.frame(alloc_strategy = c("optimumDistrict", "optimumDistrictIncid", "hrDistrict2", "hrDistrict2Incid", "optimumWat", "optimumWash", "optimumSan", "allDistrict"),
                        plot_alloc_strategy = c("case optimized", "rate optimized", "case-logistics optimized", "rate-logistics optimized", "water optimized", "watsan optimized", "sanitation optimized", "untargeted"))

  return(label_df)
}

##' generate sensitivity group labels for plotting
##' @return
label_sensitivityGroups <- function(){

  label_df <- data.frame(sensitivityLabel = c("baseline", "low_ve", "high_ve",  "low_coverage", "high_coverage", "low_campaignFreq", "high_campaignFreq", "low_vaccSupply", "high_vaccSupply", "low_lifeExpect", "high_lifeExpect", "low_indirect", "high_indirect"), 
                        plot_sensitivityLabel = c("baseline", "low vacc. efficacy", "high vacc. efficacy", "low vacc. coverage", "high vacc. coverage", "less freq. campaigns", "more freq. campaigns", "low vacc. supply", "high vacc. supply", "high pop. turnover", "low pop. turnover", "no indirect vacc. protection", "high indirect vacc. protection"),
                        sensitivityGroup = c("baseline", "vaccine efficacy", "vaccine efficacy", "campaign coverage", "campaign coverage",  "campaign frequency", "campaign frequency", "vaccine supply", "vaccine supply", "population turnover", "population turnover", "vacc. indirect protection", "vacc. indirect protection"))

  return(label_df)
}

