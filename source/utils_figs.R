
## utility functions for plotting main text figures

## district incidence in one country
plot_choro_oneCountry_lambda <- function(cntry = "COD",
                                  out_wd = "manuscripts/GAVI Impact Estimation/country_outputs_ms",
                                  fig_out_wd = "manuscripts/GAVI Impact Estimation/figures_astmh",
                                  year = 2018,
                                  w = 4, h = 4) {

  cntry_shp <- get_country_sublevels(cntry, 2)
  shp <- cntry_shp %>%
    dplyr::mutate(district = paste(NAME_0, NAME_1, NAME_2, sep = "_"))
  pop_per_district <- readRDS(paste0(out_wd,"/",cntry,"_",year,"_pop_per_district.rds"))
  cases_per_district <- readRDS(paste0(out_wd,"/",cntry,"_",year,"_cases_per_district.rds"))
  inc_per_district <- apply(cases_per_district/pop_per_district, 2, median)
  incDat <- data.frame(district = as.character(names(inc_per_district)), incid = unname(inc_per_district), log_incid = log10(unname(inc_per_district)))
  fullDat <- left_join(shp, incDat, by = c("district"))
  print(summary(pop_per_district[1,]))
  print(summary(fullDat))

  plt <- ggplot() +
    geom_sf(data = shp, colour = "grey40", fill = "grey10") +
    geom_sf(data = fullDat, aes(fill = log_incid), colour = "black") +
    scale_fill_gradientn(colours = brewer.pal(9,"Reds"), limits = log10(c(1E-6, 3E-3)), breaks = log10(c(1E-6, 1E-5, 1E-4, 1E-3)), labels = c("1M", "100K", "10K", "1K"), guide = guide_colorbar(title = expression("Incidence, 1 per")), na.value = "grey60") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.text = element_text(angle=45, vjust=1, hjust=1), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_line(colour = 'transparent'))

  ggsave(paste0(fig_out_wd, "/choro_", cntry, "_", year, "_incid_per_district.pdf"), plt, width = w, height = h, units = "in")

  return(fullDat)
}


## district targets in one country
plot_choro_oneCountry_targets <- function(cntry = "COD",
                                  alloc_strategy = "optimumDistrict",
                                  fig_out_wd = "figures/",
                                  years = 2018:2020,
                                  w = 7, h = 3,
                                  dtargets_summ) {

  cntry_shp <- get_country_sublevels(cntry, 2)
  shp <- cntry_shp %>%
    dplyr::mutate(district = paste(NAME_0, NAME_1, NAME_2, sep = "_"))

  districtSumm <- dtargets_summ %>% 
    dplyr::filter(cntry_code == cntry) %>%
    dplyr::filter(year %in% years) %>%
    dplyr::filter(!is.na(vacc_used)) %>%
    dplyr::mutate(log_vacc_used = log(vacc_used)) %>%
    dplyr::select(scenario, year, district, vacc_used, log_vacc_used) %>%
    tidyr::complete(scenario, year, district) %>%
    dplyr::filter(grepl(alloc_strategy, scenario)) %>%
    dplyr::select(-scenario)
    
  fullDat <- full_join(shp, districtSumm, sort = FALSE, by = c("district")) %>%
    dplyr::filter(!is.na(year)) %>%
    dplyr::mutate(year=factor(year, levels = years))

  plt <- ggplot() +
    geom_sf(data = shp, colour = "grey40", fill = "grey10") +
    geom_sf(data = fullDat, aes(fill = vacc_used), colour = "grey40") +
    scale_fill_gradient_tableau(palette = "Classic Blue", limits = c(2E6, 10E6), breaks = c(2E6, 4E6, 6E6, 8E6, 10E6), labels = c("2", "4", "6", "8", "10"), guide = guide_colorbar(title = expression("Doses (millions)")), na.value = "grey80") +
    facet_wrap(~year, nrow = 1) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.text = element_text(angle=45, vjust=1, hjust=1), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_line(colour = 'transparent'), strip.text = element_blank())

  ggsave(paste0(fig_out_wd, "/choro_", cntry, "_", alloc_strategy, "_dosesAlloc_", years[1], years[length(years)], ".pdf"), plt, width = w, height = h, units = "in")

  return(plt)
}


## fvp at admin 2 across all Africa
plot_choro_admin2Africa_fvp <- function(dtargets_df, allocation_strategy, 
                                        year = "allyrs", 
                                        w = 5, h = 6,
                                        fig_out_wd = "figures/",
                                        africa_shp = "Layers/Africa_SHP"){

  afr_shp <- readOGR(africa_shp, layer = "Africa")
  afr_shp2 <- st_as_sf(afr_shp)
  st_crs(afr_shp2) <- 4326

  ## sum fvps in targets across all years or choose data from a single year
  if(year == "allyrs"){
    clean_df <- dtargets_df %>%
      dplyr::filter(grepl(allocation_strategy, scenario)) %>%
      group_by(scenario, district) %>%
      dplyr::summarise(cntry_code = first(cntry_code),
                      fvp = sum(fvp, na.rm = TRUE)) %>%
      ungroup %>%
      dplyr::filter(fvp > 0)
  } else{
    clean_df <- dtargets_df %>% 
      dplyr::filter(year == year & grepl(allocation_strategy, scenario)) %>%
      dplyr::select(scenario, district, cntry_code, fvp) %>%
      dplyr::filter(fvp > 0)
  }
  cntry_ls <- unique(clean_df$cntry_code)
  
  ## for each country with any districts with non-zero fvp, compile admin2 shapefiles
  all_dshps_ls <- lapply(1:length(cntry_ls), function(i){
    cntry_shp <- get_country_sublevels(cntry_ls[i], 2)
    clean_cntry_dtargets <- clean_df %>% dplyr::filter(cntry_code == cntry_ls[i])
    shp <- cntry_shp %>%
      dplyr::mutate(district = paste(NAME_0, NAME_1, NAME_2, sep = "_")) %>%
      full_join(clean_cntry_dtargets, by = c("district"))
    return(shp)
  })
  all_dshps <- do.call(what = sf:::rbind.sf, args = all_dshps_ls)

  ## compile admin0 shapefiles
  all_shps_ls <- lapply(1:length(cntry_ls), function(i){
    cntry_shp <- get_country_shapefile(cntry_ls[i])
    clean_cntry_targets <- clean_df %>% dplyr::filter(cntry_code == cntry_ls[i])
    shp <- cntry_shp %>%
      dplyr::rename(cntry_code = ISO) %>%
      full_join(clean_cntry_targets, by = c("cntry_code"))
    return(shp)
  })
  all_shps <- do.call(what = sf:::rbind.sf, args = all_shps_ls)
  
  plt <- ggplot() +
    geom_sf(data = afr_shp2, colour = "grey60", fill = "grey90") +
    geom_sf(data = all_shps, colour = "grey60", fill = "white") +
    geom_sf(data = all_dshps, aes(fill = fvp), colour = "white", lwd = 0) +
    geom_sf(data = all_shps, colour = "grey60", fill = "transparent") +
    # coord_sf(crs = "+proj=longlat +datum=WGS84 +no_defs")
    scale_fill_viridis_c(option = "viridis", limits = c(3000, 2.7E7), breaks = c(5E3, 5E4, 5E5, 5E6), labels = c("5K", "50K", "500K", "5M"), trans = "log10", guide = guide_colorbar(title = expression("Fully Vaccinated Persons")), na.value = "white") + ##  limits = c(3000, 2.7E7) hrD2Incid & optimumDistrictIncid 
    theme_minimal() +
    theme(legend.position = "bottom", legend.text = element_text(angle=45, vjust=1, hjust=1), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_line(colour = 'transparent'))

  ggsave(paste0(fig_out_wd, "/choro_admin2Africa_", allocation_strategy, "_", year, "_fvp.pdf"), plt, width = w, height = h, units = "in")

  return(all_dshps)
}


#' @name plot_ts_cum_healthOutcomes_byScenario
#' @title plot_ts_cum_healthOutcomes_byScenario
#' @description Produce and save cumulative time series plots for heatlh outcomes across 2018-2030
#' @param gathered_df
#' @param healthOutcome
#' @param scenarioLs
#' @param w
#' @param h
#' @param fig_out_wd
#' @return 
plot_ts_cum_healthOutcomes_byScenario <- function(gathered_df, 
                                              healthOutcome, 
                                              scenarioLs, 
                                              w = 7, h = 5.5,
                                              fig_out_wd = "figures/"){
  
  pltLabels <- label_healthOutcomes() %>%
    dplyr::filter(measure == healthOutcome)
  pltLabels_allocStrategy <- label_allocStrategy() %>%
    dplyr::filter(alloc_strategy %in% unique(gathered_df$alloc_strategy))

  plot_df <- gathered_df %>%
    dplyr::filter(scenario %in% scenarioLs) %>%
    dplyr::filter(measure == healthOutcome) %>%
    dplyr::filter(year < 2031) %>%
    dplyr::mutate(year = factor(year, levels = sort(unique(gathered_df$year)))) %>%
    dplyr::mutate(alloc_strategy = factor(alloc_strategy, levels = pltLabels_allocStrategy$alloc_strategy, labels = pltLabels_allocStrategy$plot_alloc_strategy))

  plt <- ggplot(plot_df, aes(x = year, y = median, group = alloc_strategy)) +
    geom_line(aes(colour = alloc_strategy), size = 1.25) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = alloc_strategy), alpha = 0.3) + 
    scale_colour_tableau(name = "Allocation\nStrategy") +
    scale_fill_tableau(name = "Allocation\nStrategy") +
    theme_bw() + 
    theme(text = element_text(size = 14), legend.text = element_text(size = 11), legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) +
    scale_y_continuous(paste("Cumulative", pltLabels$plotMeasure), limits = c(0,NA)) +
    guides(colour=guide_legend(nrow=2))

  ggsave(paste0(fig_out_wd, "/ts_cum_", healthOutcome, "_", paste(unlist(strsplit(scenarioLs[1], "_"))[2:length(unlist(strsplit(scenarioLs[1], "_")))], collapse="_"), ".pdf"), plt, width = w, height = h, units = "in")

  return(plot_df)
}


#' @name plot_ts_relativeCases_byScenario
#' @title plot_ts_relativeCases_byScenario
#' @description Produce and save time series plot for relative cases in each vaccination deployment strategy for each year
#' @param gathered_df
#' @param scenarioLs
#' @param w
#' @param h
#' @param fig_out_wd
#' @return 
plot_ts_relativeCases_byScenario <- function(gathered_df, 
                                              scenarioLs, 
                                              w = 3.5, h = 2.5,
                                              fig_out_wd = "figures/"){
  
  pltLabels <- label_healthOutcomes() %>%
    dplyr::filter(measure %in% c("cases", "cases_cf"))
  pltLabels_allocStrategy <- label_allocStrategy() %>%
    dplyr::filter(alloc_strategy %in% unique(gathered_df$alloc_strategy))

  int_df <- gathered_df %>%
    dplyr::filter(scenario %in% scenarioLs) %>%
    dplyr::filter(measure %in% c("cases", "cases_cf")) %>%
    dplyr::mutate(alloc_strategy = factor(alloc_strategy, levels = pltLabels_allocStrategy$alloc_strategy, labels = pltLabels_allocStrategy$plot_alloc_strategy)) %>%
    dplyr::filter(year < 2031) %>%
    dplyr::mutate(year = factor(year, levels = sort(unique(gathered_df$year)))) 

  med_df <- int_df %>%
    dplyr::select(year, measure, alloc_strategy, median) %>%
    tidyr::spread(measure, median) %>%
    dplyr::mutate(med_ratio = cases/cases_cf*100) %>%
    dplyr::rename(med_cases = cases, med_cases_cf = cases_cf)
  low_df <- int_df %>%
    dplyr::select(year, measure, alloc_strategy, ci_lower) %>%
    tidyr::spread(measure, ci_lower) %>%
    dplyr::mutate(low_ratio = cases/cases_cf*100) %>%
    dplyr::rename(low_cases = cases, low_cases_cf = cases_cf)
  up_df <- int_df %>%
    dplyr::select(year, measure, alloc_strategy, ci_upper) %>%
    tidyr::spread(measure, ci_upper) %>%
    dplyr::mutate(up_ratio = cases/cases_cf*100) %>%
    dplyr::rename(up_cases = cases, up_cases_cf = cases_cf)
  plot_df <- full_join(med_df, low_df, by = c("year", "alloc_strategy")) %>%
    full_join(up_df, by = c("year", "alloc_strategy")) 

    if(cumulative == TRUE){
      plt <- ggplot(plot_df, aes(x = year, y = med_ratio, group = alloc_strategy)) +
            geom_line(aes(colour = alloc_strategy), size = 1.25) +
            geom_ribbon(aes(ymin = low_ratio, ymax = up_ratio, fill = alloc_strategy), alpha = 0.3) + 
            scale_colour_tableau(name = "Allocation\nStrategy") +
            scale_fill_tableau(name = "Allocation\nStrategy") +
            theme_bw() + 
            theme(text = element_text(size = 11), legend.text = element_text(size = 14), legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) +
            scale_y_continuous(paste0("Cumulative cases relative\nto no vaccination (%)"), limits = c(0,NA)) +
            # guides(colour=guide_legend(nrow=2))
            guides(fill = "none", colour = "none")

      ggsave(paste0(fig_out_wd, "/ts_relcum_perc_cases_", paste(unlist(strsplit(scenarioLs[1], "_"))[2:length(unlist(strsplit(scenarioLs[1], "_")))], collapse="_"), ".pdf"), plt, width = w, height = h, units = "in")  
      
     } else{
      plt <- ggplot(plot_df, aes(x = year, y = med_ratio, group = alloc_strategy)) +
            geom_line(aes(colour = alloc_strategy), size = 1.25) +
            geom_ribbon(aes(ymin = low_ratio, ymax = up_ratio, fill = alloc_strategy), alpha = 0.3) + 
            scale_colour_tableau(name = "Allocation\nStrategy") +
            scale_fill_tableau(name = "Allocation\nStrategy") +
            theme_bw() + 
            theme(text = element_text(size = 11), legend.text = element_text(size = 14), legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank()) +
            scale_y_continuous(paste0("Annual cases relative to\nno vaccination (%)", sep="\n"), limits = c(0,NA)) +
            # guides(colour=guide_legend(nrow=2))
            guides(fill = "none", colour = "none")

      ggsave(paste0(fig_out_wd, "/ts_rel_perc_cases_", paste(unlist(strsplit(scenarioLs[1], "_"))[2:length(unlist(strsplit(scenarioLs[1], "_")))], collapse="_"), ".pdf"), plt, width = w, height = h, units = "in")
      }

  return(plot_df)
}


#' @name plot_tornado_costper_oneway
#' @title plot_tornado_costper_oneway
#' @description Produce tornado plot for all sensitivity analyses
#' @param cum_ssa
#' @param scenarioLs
#' @param sensLs
#' @param meas
#' @param alloc_strategy_bl
#' @param fig_out_wd
#' @param w
#' @param h
#' @return
plot_tornado_costper_oneway <- function(cum_ssa, scenarioLs, sensLs, meas = "dalysAverted", alloc_strategy_bl = "hrDistrict2Incid", fig_out_wd, w = 5, h = 3){

  ## total median cost per FVP, discounted ##
  tcost_med <- get_vacc_cost_data(sensitivity="base")[["tot_per_FVP"]]
  dtcost_med <- get_discounted_cost_df(tcost_med)

  sensLs_blsupply <- sensLs[which((sensLs %in% c("baseline", "low_ve", "high_ve", "low_campaignFreq", "high_campaignFreq", "low_coverage", "high_coverage", "low_lifeExpect", "high_lifeExpect", "low_indirect", "high_indirect")))]
  fvp_all_df <- map_dfr(sensLs_blsupply, function(sens){
    fvp_df <- read_csv("source_inputs/vaccine_supply_baseline.csv") %>%
      dplyr::mutate(fvp = vaccines/2, sensitivityLabel = sens) %>%
      dplyr::select(sensitivityLabel, year, fvp)
    return(fvp_df)
    })
  
  fvp_hi_df <- read_csv("source_inputs/vaccine_supply_high.csv") %>%
    dplyr::mutate(fvp = vaccines/2, sensitivityLabel = "high_vaccSupply") %>%
    dplyr::select(sensitivityLabel, year, fvp) 
  fvp_lo_df <- read_csv("source_inputs/vaccine_supply_constant.csv") %>%
    dplyr::mutate(fvp = vaccines/2, sensitivityLabel = "low_vaccSupply") %>%
    dplyr::select(sensitivityLabel, year, fvp)

  fvp_df <- bind_rows(fvp_all_df, fvp_hi_df, fvp_lo_df)

  ## total median costs by sub-Saharan Africa ##
  fvps_ssa <- fvp_df %>%
    left_join(dtcost_med, by = c("year")) %>% 
    dplyr::mutate(dtcost_mn = fvp*cost) %>%
    group_by(sensitivityLabel) %>%
    summarise(dtcost_mn = sum(dtcost_mn))

  measAverted_ssa <- cum_ssa %>%
    dplyr::filter(year == 2030 & measure == meas & sensitivityLabel %in% sensLs) %>%
    full_join(fvps_ssa, by = c("sensitivityLabel")) %>%
    dplyr::select(-year) %>%
    dplyr::mutate(costper_med = ifelse(mean<0, NA, dtcost_mn/mean))

  baseline_val <- measAverted_ssa %>% 
    dplyr::filter(sensitivityLabel == "baseline" & alloc_strategy == alloc_strategy_bl) %>% 
    dplyr::select(costper_mn) %>% unlist %>% unname

  ## summary stats for campaign frequency
  cfreqSens <- measAverted_ssa %>% 
    dplyr::filter(sensitivityLabel %in% c("baseline", "high_campaignFreq", "low_campaignFreq") & alloc_strategy == alloc_strategy_bl) %>%
    dplyr::mutate(sensitivityGroup = "campaign frequency") 
  ## summary stats for coverage
  covSens <- measAverted_ssa %>% 
    dplyr::filter(sensitivityLabel %in% c("baseline", "high_coverage", "low_coverage") & alloc_strategy == alloc_strategy_bl) %>%
    dplyr::mutate(sensitivityGroup = "campaign coverage") 
  ## summary stats for ve
  veSens <- measAverted_ssa %>% 
    dplyr::filter(sensitivityLabel %in% c("baseline", "high_ve", "low_ve") & alloc_strategy == alloc_strategy_bl) %>%
    dplyr::mutate(sensitivityGroup = "vaccine efficacy") 
  ## summary stats for targeting
  targetSens <- measAverted_ssa %>%
    dplyr::filter(sensitivityLabel == "baseline" & scenario %in% scenarioLs) %>%
    dplyr::mutate(sensitivityGroup = "deployment strategy") 
  ## summary stats for vaccine supply
  vsupplySens <- measAverted_ssa %>%
    dplyr::filter(sensitivityLabel %in% c("baseline", "high_vaccSupply", "low_vaccSupply") & alloc_strategy == alloc_strategy_bl) %>%
    dplyr::mutate(sensitivityGroup = "vaccine supply") 
  ## summary stats for population turnover rate
  lfexpectSens <- measAverted_ssa %>%
    dplyr::filter(sensitivityLabel %in% c("baseline", "high_lifeExpect", "low_lifeExpect") & alloc_strategy == alloc_strategy_bl) %>%
    dplyr::mutate(sensitivityGroup = "population turnover") 
  ## summary stats for vaccine indirect protection
  indirectSens <- measAverted_ssa %>%
    dplyr::filter(sensitivityLabel %in% c("baseline", "high_indirect", "low_indirect") & alloc_strategy == alloc_strategy_bl) %>%
    dplyr::mutate(sensitivityGroup = "vacc. indirect protection")

  plotdat <- bind_rows(cfreqSens, covSens, veSens, targetSens, vsupplySens, lfexpectSens, indirectSens) %>%
    group_by(sensitivityGroup) %>%
    summarise(costper_min = min(costper_mn), costper_max = max(costper_mn)) %>%
    ungroup %>%
    dplyr::mutate(costper_min_sc = costper_min/baseline_val, costper_max_sc = costper_max/baseline_val) %>%
    dplyr::mutate(range_diff = costper_max-costper_min) %>%
    dplyr::arrange(range_diff) %>%
    dplyr::mutate(sensitivityGroup = factor(sensitivityGroup, levels = sensitivityGroup))
  
  rawplotdat <- bind_rows(cfreqSens, covSens, veSens, targetSens, vsupplySens, lfexpectSens, indirectSens) %>%
    group_by(sensitivityGroup) %>%
    dplyr::mutate(minval = ifelse(min(costper_mn)==costper_mn, TRUE, FALSE),
                  maxval = ifelse(max(costper_mn)==costper_mn, TRUE, FALSE)) %>%
    ungroup
  untargeted <- rawplotdat %>% dplyr::filter(sensitivityGroup == "deployment strategy" & alloc_strategy == "allDistrict") %>% dplyr::select(costper_mn) %>% unlist %>% unname
  rateOptim <- rawplotdat %>% dplyr::filter(sensitivityGroup == "deployment strategy" & alloc_strategy == "optimumDistrictIncid") %>% dplyr::select(costper_mn) %>% unlist %>% unname
  watOptim <- rawplotdat %>% dplyr::filter(sensitivityGroup == "deployment strategy" & alloc_strategy == "optimumWat") %>% dplyr::select(costper_mn) %>% unlist %>% unname
  sanOptim <- rawplotdat %>% dplyr::filter(sensitivityGroup == "deployment strategy" & alloc_strategy == "optimumSan") %>% dplyr::select(costper_mn) %>% unlist %>% unname

  plt <- ggplot(plotdat, aes(x = sensitivityGroup)) +
    geom_linerange(aes(ymin = costper_min, ymax = costper_max), colour = "#008080", size = 5) +
    annotate("segment", y = untargeted, yend = untargeted, x = 6.78, xend = 7.22, colour = "grey60", size=.2) + ## untargeted
    annotate("segment", y = rateOptim, yend = rateOptim, x = 6.78, xend = 7.22, colour = "grey60", size=.2) + ## rate optimized
    annotate("segment", y = watOptim, yend = watOptim, x = 6.78, xend = 7.22, colour = "grey60", size=.2) + ## water optimized
    annotate("segment", y = sanOptim, yend = sanOptim, x = 6.78, xend = 7.22, colour = "grey60", size=.2) + ## sanitation optimized
    geom_hline(aes(yintercept = baseline_val), colour = "#800000") +
    coord_flip() +
    theme_bw() +
    theme(text = element_text(size = 12), axis.title.y = element_blank(), panel.grid = element_blank()) +
    scale_y_continuous("Cost per DALY averted (2017 USD)", expand = expand_scale(mult = .3), trans = "log10")

  ggsave(paste0(fig_out_wd, "/tornado_costper_", meas, "_sensitivity_v2.pdf"), plt, width = w, height = h, units = "in")
  return(measAverted_ssa)

}

