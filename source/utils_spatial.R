
## these functions represent raster and spatial data utilities employed by the model
## In the near future, several of the shapefile-related functions will be replaced with more stable and tested versions from an R package soon-to-be deployed on CRAN (last updated September 2019)

library(R6)
library(Matrix)
library(dplyr)
library(lubridate)
library(lwgeom)
library(raster)
library(reshape2)
library(sf)
library(tidyr)
library(units)
library(sp)
library(rgdal)


#' @export R6raster
#' @title R6raster
#' @description A wrapper for the various Raster* classes from the raster package.  This is just done to prevent pass by copy.
#' @docType class
#' @importFrom R6 R6Class
R6raster <- R6Class(
  classname="R6raster",
  private = list(
    #' @importFrom raster raster
    .raster = raster()
  ),
  public = list(
    initialize = function(raster_object){
      self$raster = raster_object
    },
    load_from_file = function(filename){
      if(length(filename) == 1){
        #' @importFrom raster brick
        private$.raster = raster::brick(filename)
      } else {
        #' @importFrom raster stack
        private$.raster = raster::stack(filename)
      }
    }
  ),
  active = list(
    raster = function(value){
      if(missing(value)){
        return(private$.raster)
      }
      if(any(c('RasterLayer','RasterBrick','RasterStack') %in% class(value))){
        private$.raster = value
        return()
      }
      if(mode(class(value)) == 'character'){
        self$load_from_file(value)
        return()
      }
      stop("raster must be a RasterLayer, RasterStack, or RasterBrick object")
    }
  )
)


#' @export
#' @name lookup_WorldPop_region
#' @title lookup_WorldPop_region
#' @description Determine which WorldPop raster to look at for population of a country
#' @param location The name of the country
#' @param verbose Whether to print detailed error messages
lookup_WorldPop_region <- function(location, verbose=TRUE){
  if('character' %in% class(location)){
    fname <- system.file("extdata","country_worldpop_regions.csv",package = "taxdat")
    if(nchar(fname) == 0){
      if(verbose){
        warning("Could not load country_worldpop_regions.csv from the package directory.  Try reinstalling the package or contacting the maintainer.")
      }
      return(NA)
    }
    worldpop_regions = tryCatch(
      read.csv(
        fname,
        sep=',',
        header=TRUE,
        stringsAsFactors=FALSE,
        na.strings = "",
        colClasses = 'character',
        quote = "\"",
        check.names = FALSE
      ),
      warning = function(w){
        if(length(w$message > 0)){
          if(grepl(pattern="ncomplete final line",x=w$message)){
            return(suppressWarnings(
              read.csv(
                fname,
                sep=',',
                header=TRUE,
                stringsAsFactors=FALSE,
                na.strings = "",
                colClasses = 'character',
                quote = "\"",
                check.names = FALSE
              )
            ))
          }
        }
        if(verbose){
          warning(paste(fname,":",w$message),immediate.=TRUE)
        }
        return(w)
      },
      error = function(e){
        if(verbose){
          warning(paste(fname,":",e),immediate.=TRUE)
        }
        return(e)
      }
    )
    if(location %in% worldpop_regions[,'country']){
      return(worldpop_regions[location == worldpop_regions[,'country'],'region_code'])
    }
    if(location %in% worldpop_regions[,'country_code']){
      return(worldpop_regions[location == worldpop_regions[,'country_code'],'region_code'])
    }
  
    if(verbose){
      warning(paste("Could not find WorldPop data for country",location))
      return(NA)
    }
  } else if('sf' %in% class(location)){
    all_Worldpop_regions <- c("AMR","AFR","OCE","ASI")
    intersects = c()
    for(region in all_Worldpop_regions){
      files <- list.files(paste("data",region,sep='/'))
      # files <- files[grepl('adj',files)]
      files <- files[grepl('tif$',files)]
      files <- files[1]
      raster_layer <- raster(paste("data",region,files,sep='/'))
      # #' @importFrom sf st_bbox
      # raster_shapefile<- st_bbox(raster_layer)
      #' @importFrom sf st_intersects
      intersects[region] = st_intersects(
        #' @importFrom sf st_as_sfc
        #' @importFrom sf st_bbox
        st_as_sfc(st_bbox(raster_layer)),
        #' @importFrom sf st_as_sfc
        #' @importFrom sf st_bbox
        st_as_sfc(st_bbox(location$geometry)),
        sparse = FALSE
      )
    }
    if(sum(intersects) == 1){
      return(names(intersects)[intersects])
    } else {
      stop("Not yet written")
    }
  }
}


#' @name load_population_estimates
#' @title load_population_estimates
#' @description Find population estimates for a year and overlay them on a raster.  This will take much longer on rasters without a calculable resolution
#' @param region The WorldPop region for which underlying raster to use. "AFR" 
#' @param raster A R6raster object used to determine the extent and resolution
#' @param year Which year to try to estimate
#' @param layers_dir Layers directory of the taxonomy folder
#' @return A RasterLayer object of the same resolution and extent as \code{raster}
load_population_estimates <- function(region, raster_layer, year, layers_dir = "Layers"){
  if(!("R6raster" %in% class(raster_layer))){ stop("load_population_estimates raster_layer must be of class R6raster")}
  allowed_years = c(2000,2005,2010,2015,2020)
  
  if(year > max(allowed_years)){
    warning(paste("Data used to calculate population only extends to ",max(allowed_years)))
  }
  
  if(year < min(allowed_years)){
    warning(paste("Data used to calculate population starts from ",min(allowed_years)))
  }
  ## Get all of the filenames for population rasters
  fnames <- list(
    AFR = paste("data/AFR_PPP",
            allowed_years,
            "adj_v2.tif",
            sep='_')
  )
  
  ## if it is, set it to pop_file
  fname <- fnames[[region]]
  if(!all(file.exists(fname))){
    warning("Could not load population estimates for the following years:",paste(allowed_years[!file.exists(fname)]))
    allowed_years = allowed_years[file.exists(fname)]
    fname = fname[file.exists(fname)]
    if(length(fname) < 2){stop("Not enough population estimates found")}
  }
  #' @importFrom raster stack
  pop_raster <- R6raster$new(fname)
  #' @importFrom raster intersect
  if(is.null(raster::intersect(
    extent(pop_raster$raster),
    extent(raster_layer$raster)
  ))){
    message("Copying raster because all NA")
    pop_raster <- raster_layer$clone(TRUE)
    pop_raster$raster[] <- NA
  }
  
  ##This should probably be replaced with fun=sum
  #' @importFrom raster aggregate
  pop_raster$raster <- aggregate(pop_raster$raster,resol[1])*resol[1]*resol[2]
  if(all(is.na(pop_raster$raster[]))){
    stop("The region does not contain the raster.")
  }
  #' @importFrom raster values values<-
  
  ## should do this: set pops to NA for cells with NA in the every layer of the input raster
  ## preferable to pass r6 object to updated taxdat functions
  values(pop_raster$raster)[aggregate_raster_xlayers(raster_layer,function(x){all(is.na(x))})$raster[],] <- NA
  
  #' @importFrom reshape2 melt
  #' @importFrom raster values
  tmp = t(as.matrix(values(pop_raster$raster))) ## mx of 5 years of worldpop data, years as columns
  non_na_indices = apply(tmp,2,function(x){!all(is.na(x))})
  tmp = tmp[,non_na_indices] ## rm gridcells that are masked and where pop is NA
  models <- apply(tmp,2,function(x){
    lm(log(1+y)~x,data.frame(y=x,x=allowed_years)) ## estimate and project population
  })
  
  ## create new raster with one layer
  pop_raster <- aggregate_raster_xlayers(pop_raster,function(x){x[[1]]})
  #' @importFrom raster values values<-
  values(pop_raster$raster)[non_na_indices] = exp(sapply(models,predict,newdata=data.frame(x=year)))-1
  return(pop_raster)
}


#' @name extract_country_from_raster
#' @title extract_country_from_raster
#' @description Extract a R6raster data only for locations in a shapefile
#' @param raster_layer The R6raster to extract from.
#' @param shapefile The shapefile to use for bounding
#' @return A Raster* of the same type as \code{raster_layer}
#' @param partial_cover whether or not to do partial matching
#' @param method Whether to use the raster method or the previous hard coded one
#' @param trim_small Remove values with less than this fraction of the median overlap
extract_country_from_raster <- function(raster_layer, shapefile, partial_cover = FALSE, method='raster', trim_small=0){
  if(!("R6raster" %in% class(raster_layer))){
    raster_layer = R6raster$new(raster_layer)
  }
  if("sf" %in% class(shapefile)){
    shapefile = as_Spatial(shapefile$geometry)
  }
  if("sfc" %in% class(shapefile)){
    #' @importFrom sf as_Spatial
    shapefile = as_Spatial(shapefile)
  }
  #' @importFrom raster extent
  lhs = extent(raster_layer$raster)
  #' @importFrom raster extent
  ## This should maybe be imported from sf somehow?
  #' @importFrom raster extent
  rhs = extent(shapefile)
  if(!compare_extents(shapefile,raster_layer)){
    #' @importFrom raster crop
    message("Copying raster because extents are not equal.")
    raster_layer <- R6raster$new(crop(
      raster_layer$raster,
      #' @importFrom raster extent
      extent(shapefile),
      snap='out'
    ))
  }
  if(partial_cover){
    if(method == 'raster'){
      if(trim_small <= 0.01){
        stop("For method raster, trim_small must be at least .01")
      }
      #' @importFrom raster rasterize
      tmp_raster <- rasterize(shapefile,raster_layer$raster,getCover=partial_cover)/100
      tmp_raster[tmp_raster == 0] <- NA
      #' @importFrom raster values values<-
      values(tmp_raster)[values(tmp_raster) < quantile(values(tmp_raster),.5,na.rm=T)*trim_small] <- NA
      return(R6raster$new(tmp_raster * raster_layer$raster))
    } else {
      if(trim_small != 0){
        stop("Not written")
      }
      na_only = FALSE
      #' @importFrom raster values values<-
      if(all(is.na(unique(values(raster_layer$raster))))){
        message("Copying raster because it is all NA")
        raster_layer = raster_layer$clone(TRUE)
        na_only = TRUE
        ## FIX ME
        values(raster_layer$raster) <- 0
      }
      #' @importFrom sf st_union
      #' @importFrom sf st_as_sf
      shapefile = st_union(st_as_sf(shapefile))
      #' @importFrom raster raster
      #' @importFrom sf st_as_sf
      all_polys = st_as_sf(rasterToPolygons(raster_layer$raster,dissolve = FALSE))
      #' @importFrom sf st_crs
      st_crs(all_polys) <- st_crs(shapefile)
      #' @importFrom sf st_intersects
      intersects = sapply(st_intersects(all_polys,shapefile),length) > 0
      if(na_only){
        intersects = is.na(intersects)
      }
      #' @importFrom sf st_intersection
      cropped_polys = st_intersection(all_polys[intersects,],shapefile)
      #' @importFrom sf st_area
      adjustment_factors = st_area(cropped_polys)/st_area(all_polys[intersects,])
      if(class(raster_layer$raster) == "RasterBrick"){
        indices_1 <- which(apply(raster_layer$raster@data@values,1,function(x){any(!is.na(x))}))
        indices_2 = indices_1[intersects]
        indices_3 = indices_1[!(indices_1 %in% indices_2)]
        message("Copying raster to make return object")
        raster_layer = raster_layer$clone(TRUE)
        raster_layer$raster@data@values[indices_2,] = raster_layer$raster@data@values[indices_2,] * adjustment_factors
        raster_layer$raster@data@values[indices_3,] = NA
      } else if(class(raster_layer$raster) == "RasterStack"){
        message("Copying raster to make return object")
        raster_layer = raster_layer$clone(TRUE)
        #' @importFrom raster values
        indices_1 <- which(apply(values(raster_layer$raster),1,function(x){any(!is.na(x))}))
        indices_2 = indices_1[intersects]
        indices_3 = indices_1[!(indices_1 %in% indices_2)]
        #' @importFrom raster values values<-
        values(raster_layer$raster)[indices_2,] = values(raster_layer$raster)[indices_2,] * adjustment_factors
        #' @importFrom raster values  values<-
        values(raster_layer$raster)[indices_3,] = NA
      } else {
        message("Copying raster to make return object")
        raster_layer = raster_layer$clone(TRUE)
        indices_1 <- which(!is.na(raster_layer$raster@data@values))
        indices_2 = indices_1[intersects]
        indices_3 = indices_1[!(indices_1 %in% indices_2)]
        raster_layer$raster@data@values[indices_2] = raster_layer$raster@data@values[indices_2] * adjustment_factors
        raster_layer$raster@data@values[indices_3] = NA
      }
    }
  } else {
    message("Copying raster to make return object")
    raster_layer = raster_layer$clone(TRUE)
    #' @importFrom raster mask
    raster_layer$raster <- mask(
      raster_layer$raster,
      shapefile
    )
  }
  
  return(raster_layer)
  # #' @importFrom raster extract
  # extract(raster_layer,shapefiles)
}


#' @name aggregate_raster_xlayers
#' @title aggregate_raster_xlayers
#' @description A mask of stackApply that always assumes indices is all 1
#' @param x Raster* object
#' @param fun function that returns a single value, e.g. mean or min, and that takes a na.rm argument (or can pass through arguments via ...)
#' @param ... additional arguments as for writeRaster
#' @param sections integer As part of aggregating, break the raster into sections x sections pieces first.
#' @return A RasterLayer
aggregate_raster_xlayers <- function(x, fun, ..., method = 'apply', sections=1){
  if(!("R6raster" %in% class(x))){
    x = R6raster$new(x)
  }
  #' @importFrom raster nlayers
  if(nlayers(x$raster) <= 1){
    return(x)
  }
  rc <- R6raster$new(x$raster[[1]])
  ## Making sure people dont go overboard on big rasters.
  if(sections == 1){
    sections = max(round(dim(x$raster)/1000),1)
    if(sections != 1){warning(paste("This raster is large, breaking it up into",sections,"sections"))}
  }

  row_starts = floor(seq(1,dim(x$raster)[1],(dim(x$raster)/sections)[1]))
  row_ends = c(row_starts[-1] - 1,dim(x$raster)[1])
  row_lengths = row_ends - row_starts + 1
  col_starts = floor(seq(1,dim(x$raster)[2],(dim(x$raster)/sections)[2]))
  col_ends = c(col_starts[-1] - 1,dim(x$raster)[2])
  col_lengths = col_ends - col_starts + 1
  #' @import foreach
  foreach(row_section = 1:length(row_starts)) %do% {
    foreach(col_section = 1:length(col_starts)) %do% {
      rc$raster[row_starts[row_section]:row_ends[row_section],col_starts[col_section]:col_ends[col_section]] <-
        #' @importFrom raster getValuesBlock
        apply(getValuesBlock(x$raster,row_starts[row_section],row_lengths[row_section],col_starts[col_section],col_lengths[col_section]),1,fun,...)
    }
  }
  return(rc)
}


#' @name apply_to_all_sublevels
#' @title apply_to_all_sublevels
#' @description Apply a function to a raster with different answers for each subregion of a country at a particular ISO level
#' @param location_name The name of the country
#' @param raster The raster to apply the function to
#' @param ISO_level The level to get the shapefiles at
#' @param taxonomy_dir The taxonomy directory to store the shapefiles in
#' @param verbose Set to true for more warnings and messages
#' @param partial_cover Set to true to get accurate partial matches by subsampling the raster grid.
apply_to_all_sublevels <- function(location_name,
                                   raster_layer,
                                   ISO_level,
                                   fun,
                                   taxonomy_dir = 'taxonomy-verified',
                                   verbose=FALSE,
                                   fix_location_names = FALSE,
                                   partial_cover = FALSE,
                                   trim_small = 0,
                                   method = 'not raster',
                                   ...){
  
  if(!("R6raster" %in% class(raster_layer))){
    raster_layer = R6raster$new(raster_layer)
  }
  
  all_shapefiles <- get_country_sublevels(location_name,ISO_level,taxonomy_dir,verbose)
  rc <- apply(
    all_shapefiles,
    1,
    function(shp){
      shp = st_sfc(shp$geometry)
      extracted_raster = extract_country_from_raster(
        raster_layer,
        shp,
        partial_cover = partial_cover,
        trim_small = trim_small,
        method = method
      )
      aggregate_raster_xcells(extracted_raster,fun,...)
    }
  )
  
  if(fix_location_names){
    rc_names <- apply(
      sapply(
        0:ISO_level,
        function(level){
          all_shapefiles[[paste("NAME",level,sep="_")]]
        }
      ),
      1,
      function(name){
        fix_location_name(
          paste('???',paste(name,collapse='_'),sep='_'),
          taxonomy_dir = taxonomy_dir,
          verbose=verbose
        )
      }
    )
  } else {
    rc_names <- apply(
      sapply(
        0:ISO_level,
        function(level){
          all_shapefiles[[paste("NAME",level,sep="_")]]
        }
      ),
      1,
      paste,
      collapse='_'
    )
  }
  # rc_names <- all_shapefiles[[paste("NAME",1:ISO_level,sep='_')]]
  if(is.null(dim(rc))){
    names(rc) <- rc_names
  } else {
    colnames(rc) <- rc_names
  }
  return(rc)
}


#' @name aggregate_raster_xcells
#' @title aggregate_raster_xcells
#' @description Apply a function to each layer, aggregating the cells of the raster
#' @param x Raster* object
#' @param fun function that returns a single value, e.g. mean or min, and that takes a na.rm argument (or can pass through arguments via ...)
#' @param ... additional arguments as for writeRaster
#' @return A vector
aggregate_raster_xcells <- function(x,fun,...){
  if(!("R6raster" %in% class(x))){
    x = R6raster$new(x)
  }
  if(is.null(dim(getValues(x$raster)))){
    #' @importFrom raster getValues
    return(fun(getValues(x$raster),...))
  }
  #' @importFrom raster getValues
  apply(getValues(x$raster),2,fun,...)
}


#' @name compare_extents
compare_extents <- function(shapefile,raster_layer){
  #' @importFrom raster extent
  lhs = extent(shapefile)
  #' @importFrom raster extent
  if("R6raster" %in% class(raster_layer)){
    rhs = extent(raster_layer$raster)
    resol = res(raster_layer$raster)
  } else{
    rhs = extent(raster_layer)
    resol = res(raster_layer)
  }
  #' @importFrom raster res  
  return(
    (abs(lhs@xmin - rhs@xmin) < resol[1]) &
      (abs(lhs@xmax - rhs@xmax) < resol[1]) &
      (abs(lhs@ymin - rhs@ymin) < resol[2]) &
      (abs(lhs@ymax - rhs@ymax) < resol[2])
  )
}


#' @name crop_raster_to_shapefile
#' @title crop_raster_to_shapefile
#' @description Crop a raster to the extent of a shapefile
#' @param raster_layer The R6raster to extract from.
#' @param shapefile The shapefile to use for bounding
#' @return A R6raster of the same type as \code{raster_layer} with the extent determined by \code{shapefile}
crop_raster_to_shapefile <- function(raster_layer,shapefile){
  if(!("R6raster" %in% class(raster_layer))){
    raster_layer = R6raster$new(raster_layer)
  }
  if("sf" %in% class(shapefile)){
    shapefile = as_Spatial(shapefile$geometry)
  }
  if("sfc" %in% class(shapefile)){
    #' @importFrom sf as_Spatial
    shapefile = as_Spatial(shapefile)
  }
  ## This should maybe be imported from sf somehow?
  #' @importFrom raster extent
  if(!compare_extents(shapefile,raster_layer$raster)){
    #' @importFrom raster crop
    message("Copying raster because extents are not equal.")
    raster_layer <- R6raster$new(crop(
      raster_layer$raster,
      #' @importFrom raster extent
      extent(shapefile),
      snap='out'
    ))
  }
  return(raster_layer)
}


#' @name get_country_sublevels
#' @title get_country_sublevels
#' @param location_name The name of the country to get a shapefile for
#' @param ISO_level What ISO level to pull shapefiles for.  1 is the level immediately below the country level
#' @param taxonomy_dir The path to the taxonomy filesystem to store the shapefile in.
#' @param land_only boolean If TRUE, try to filter out shapefiles that are not associated with land regions.
#' @param verbose Set to TRUE for more warnings and messages.
#' @importFrom lubridate now
get_country_sublevels = function(location_name,ISO_level,taxonomy_dir = 'taxonomy-verified', land_only=FALSE,method='name',layers_dir = 'Layers',date=now(),verbose=method=='center'){
  location_name = fix_location_name(paste("???",fix_country_name(location_name),sep='_'),taxonomy_dir = taxonomy_dir,verbose=FALSE)
  ISO_A1 = strsplit(x = location_name,split="_")[[1]][2]
  ## Check for custom shapefiles.
  
  if(method == 'center'){
    location_dir <- paste(taxonomy_dir,'Location',sep='/')
    alllocs <- list.files(location_dir,pattern = location_name,ignore.case = TRUE)
    # pattern = paste('^[^_]*',paste(rep('_',ISO_level + 2),collapse='[^_]*'),'LOC.csv',sep='')
    # alllocs = alllocs[grep(pattern,alllocs)]
    allshp <- lapply(
      alllocs,
      function(loc){
        locfile <- read.csv(paste(taxonomy_dir,'Location',loc,sep='/'),stringsAsFactors = FALSE,header = FALSE)
        iso_loc_indices = grep('isISO_A2_L',locfile[,1])
        iso_loc_indices = iso_loc_indices[as.numeric(locfile[iso_loc_indices,2]) %in% c(1,2)]
        iso_level <- max(c(0,as.numeric(substr(locfile[iso_loc_indices,1],nchar(locfile[iso_loc_indices,1]),nchar(locfile[iso_loc_indices,1])))))
        if(iso_level != ISO_level){
          return(NULL)
        }
        loc_name = strsplit(loc,'_')[[1]]
        loc_name = loc_name[1:(length(loc_name)-1)]
        loc_name = paste(loc_name,collapse='_')
        return(get_shape_file(location_name = loc_name,taxonomy_dir = taxonomy_dir,layers_dir = layers_dir,verbose = verbose,method = 'center',date = date))
      }
    )
    if(any(duplicated(allshp) & (!sapply(allshp, is.null)))){
      if(verbose){
        duplicate_idx <- duplicated(allshp) & (!sapply(allshp, is.null))
        for(idx in duplicate_idx){
          matches = sapply(allshp,function(x){all.equal(allshp[[idx]],x) == TRUE})
          paste('Location',alllocs[idx],'is a duplicate of ',alllocs[matches])
        }
      }
    }
    allshp <- reduce_sf_vector(unique(allshp))
    if(verbose){
      cntry_shp <- get_shape_file(location_name)
      #' @importFrom sf st_equals
      if(!st_equals(cntry_shp$geometry,st_union(allshp),sparse = FALSE)){
        warning("Sublevels incomplete")
      }
    }
    return(allshp)
  }

  destination = paste(taxonomy_dir,'/ShapeFiles/',ISO_A1,'_adm',ISO_level,'.rds',sep='')

  if(!file.exists(paste(taxonomy_dir,"ShapeFiles",sep='/'))){
    stop("The directory ",taxonomy_dir," should exist, and have a subdirectory called ShapeFiles")
  }
  if(!file.exists(destination)){
    website = paste("http://biogeo.ucdavis.edu/data/gadm2.8/rds/",ISO_A1,"_adm",ISO_level,".rds",sep='')
    download.file(website,destination,mode='wb')
  }
  #' @importFrom sf st_as_sf
  all_shape_files = st_as_sf(readRDS(destination))
  #' @importFrom sf st_is_valid
  if(!all(suppressWarnings(st_is_valid(all_shape_files)))){
    #' @importFrom lwgeom st_make_valid
    all_shape_files <- st_make_valid(all_shape_files)
  }
  #' @importFrom sf st_is_valid
  if(!all(suppressWarnings(st_is_valid(all_shape_files)))){
    #' @importFrom lwgeom st_make_valid
    all_shape_files <- st_make_valid(all_shape_files)
  }
  if(land_only){
    all_shape_files <- all_shape_files[all_shape_files[[paste("TYPE",ISO_level,sep='_')]] != "Water body",]
  }
  return(all_shape_files)
}


#' @export
#' @name fix_country_name
#' @title fix_country_name
#' @description Change a country name to its ISO_3166_1 alpha-3 code
#' @param country_name The name of the country
fix_country_name <- function(country_name,verbose=TRUE){
  fname <- system.file("extdata","country_aliases.csv",package = "taxdat")
  if(nchar(fname) == 0){
    if(verbose){
      warning("Could not load country_aliases.csv from the package directory.  Try reinstalling the package or contacting the maintainer.")
    }
    return(NA)
  }
  country_aliases  = tryCatch(
    read.csv(
      fname,
      sep=',',
      header=TRUE,
      stringsAsFactors=FALSE,
      na.strings = "",
      colClasses = 'character',
      quote = "\"",
      check.names = FALSE
    ),
    warning = function(w){
      if(length(w$message > 0)){
        if(grepl(pattern="ncomplete final line",x=w$message)){
          return(suppressWarnings(
            read.csv(
              fname,
              sep=',',
              header=TRUE,
              stringsAsFactors=FALSE,
              na.strings = "",
              colClasses = 'character',
              quote = "\"",
              check.names = FALSE
            )
          ))
        }
      }
      if(verbose){
        warning(paste(fname,":",w$message),immediate.=TRUE)
      }
      return(w)
    },
    error = function(e){
      if(verbose){
        warning(paste(fname,":",e),immediate.=TRUE)
      }
      return(e)
    }
  )
  if(country_name %in% country_aliases[,2]){
    return(country_aliases[country_name == country_aliases[,2],1])
  } else if(country_name %in% country_aliases[,1]){
    return(country_name)
  } 
  
  if(verbose){
    warning(paste("Could not find a country with the name",country_name))
    return(country_name)
  } 
}


#' @name get_country_shapefile
#' @title get_country_shapefile
#' @param name The name of the country to get a shapefile for
#' @param taxonomy_dir The path to the taxonomy filesystem to store the shapefile in.
get_country_shapefile = function(name,taxonomy_dir = 'taxonomy-verified',verbose=FALSE){
  location_name = paste("???",fix_country_name(name),sep='_')
  if(!is.na(location_name)){
    return(get_shape_file(location_name,taxonomy_dir,verbose=FALSE))
  }
  return(NA)
}


#' @name get_shape_file
#' @title get_shape_file
#' @description Look up the GADM shapefile for a particular location
#' @param location Name of a location as a character
#' @param taxonomy_dir Name of a taxonomy directory (shapefiles will be stored here)
#' @param layers_dir Name of a layers directory (custom shapefiles will be found here)
#' @param method which method to use when pulling shapefiles.  If \code{method =='center'}, then use the central points listed in the shape files.  If \code{method=='name'} use the names.
#' @return an sf::sf with the results
#' @importFrom lubridate now
get_shape_file <- function(location_name,taxonomy_dir = 'taxonomy-verified',layers_dir = 'Layers', verbose=TRUE,method='name',date=now()){
  ##Check to see if a custom shapefile exists:
  custom_shape_file <- exists_shape_file(location_name,taxonomy_dir,verbose,date=date)
  using_custom <- custom_shape_file != FALSE
  if(verbose && ('_' == substr(location_name,nchar(location_name),nchar(location_name)))){
    warning("Location ",location_name," contains an ending underscore.")
  }
  rc_shape_files <- NULL
  if(verbose && ('_' == substr(location_name,nchar(location_name),nchar(location_name)))){
    warning("Location ",location_name," contains an ending underscore.")
  }

  location.tmp = location_name
  location.tmp = gsub('-','dashcharacter',location.tmp)
  location.tmp = gsub('|','vertcharacter',location.tmp,fixed=T)
  location.tmp = gsub('_','underscorecharacter',location.tmp)
  location.tmp = iconv(location.tmp, from="UTF-8",to="ASCII//TRANSLIT")
  location.tmp = gsub(' ','',location.tmp)
  location.tmp = gsub("[[:punct:]]",'',location.tmp)
  location.tmp = gsub('dashcharacter','-',location.tmp)
  location.tmp = gsub('vertcharacter','|',location.tmp)
  location.tmp = gsub('underscorecharacter','_',location.tmp)
  location.tmp = gsub('__','_',location.tmp)
  location.tmp = gsub('_$','',location.tmp)
  
  ## Read location files
  location_files <- list.files(paste(taxonomy_dir,'Location',sep='/'))
  original_location_files <- location_files
  location_files = gsub('-','dashcharacter',location_files)
  location_files = gsub('|','vertcharacter',location_files,fixed=T)
  location_files = gsub('_','underscorecharacter',location_files)
  location_files = iconv(location_files, from="UTF-8",to="ASCII//TRANSLIT")
  location_files = gsub(' ','',location_files)
  location_files = gsub("[[:punct:]]",'',location_files)
  location_files = gsub('dashcharacter','-',location_files)
  location_files = gsub('vertcharacter','|',location_files)
  location_files = gsub('underscorecharacter','_',location_files)
  location_files = gsub('__','_',location_files)
  location_files = gsub('_$','',location_files)
  location_files = sapply(
    strsplit(location_files,'_'),
    function(x){paste(x[-length(x)],collapse = '_')}
  )
  location_files = toupper(location_files)
  location_files = gsub('-','',location_files)
  location_files <- setNames(original_location_files,location_files)

  if(verbose && (location.tmp != location_name)){
    warning("Location ",location_name," contains illegal characters")
  }
  location_name = toupper(location.tmp)
  location_name = gsub('-','',location_name)

  location.array = strsplit(location_name,split='_')[[1]]
  who_region = location.array[[1]]
  if(length(location.array) == 1){
    # browser()
    stop("Not yet written")
  }
  all_ISO_A1 = strsplit(location.array[[2]],"|",fixed=TRUE)[[1]]
  ISO_A2 = location.array[-c(1,2)]
  ISO_A2_level = length(ISO_A2)
  all_ISO_A2 = strsplit(ISO_A2,"|",fixed=TRUE)
  if(ISO_A2_level > 3){
    if(verbose){
      message("GADM does not normally have ISO levels this deep")
    }
  }
  for(loc_idx in 1:length(all_ISO_A1)){
    ISO_A1 = all_ISO_A1[loc_idx]
    ISO_A2 = sapply(all_ISO_A2,function(x){x[loc_idx]})

    if(using_custom){
      #' @importFrom sf st_read
      all_shape_files = st_read(paste(layers_dir,'shapefiles',custom_shape_file,sep='/'),quiet=TRUE)
      method = 'center'
    } else {
      destination = paste(taxonomy_dir,'/ShapeFiles/',ISO_A1,'_adm',ISO_A2_level,'.rds',sep='')
      if(!file.exists(destination)){
        website = paste("http://biogeo.ucdavis.edu/data/gadm2.8/rds/",ISO_A1,"_adm",ISO_A2_level,".rds",sep='')
        download.file(website,destination,mode='wb')
      }
      #' @importFrom sf st_as_sf
      all_shape_files = st_as_sf(readRDS(destination))
    }
    #' @importFrom sf st_is_valid
    if(!all(suppressWarnings(st_is_valid(all_shape_files)))){
      #' @importFrom lwgeom st_make_valid
      all_shape_files <- st_make_valid(all_shape_files)
    }
    if(method == 'center'){
      ## Lookup location file
      this_location_file = paste(
        taxonomy_dir,
        "Location",
        location_files[gsub('_$','',paste(who_region,ISO_A1,paste(ISO_A2,collapse='_'),sep='_'))],
        sep='/'
      )
      if(sapply(strsplit(this_location_file,'/'),function(x){x[length(x)]}) == "NA"){
        if(verbose){
          warning(paste("Could not find location file for location",location_name))
        }
        sublocation_files <- location_files[startsWith(names(location_files),paste(who_region,ISO_A1,paste(ISO_A2,collapse='_'),sep='_'))]
        if(length(sublocation_files) > 0){
          if(verbose){
            warning(paste("Using", length(sublocation_files),"sub-locations as a proxy for",location_name))
          }
          loc_data <- lapply(paste(taxonomy_dir,'Location',sublocation_files,sep='/'),read_location_csv)
          #' @importFrom sf st_multipoint
          central_points <- st_multipoint(matrix(c(sapply(loc_data,function(x){x$cent_long}),sapply(loc_data,function(x){x$cent_lat})),length(loc_data),2))
          #' @importFrom sf st_contains
          this_shape_file <- all_shape_files[st_contains(all_shape_files,central_points,sparse=F),]
        } else {
          this_shape_file <- all_shape_files[0,]
        }
      } else {
        loc_data <- read_location_csv(this_location_file)
        #' @importFrom sf st_point
        loc_point <- st_point(c(loc_data$cent_long,loc_data$cent_lat))
        #' @importFrom sf st_intersects
        this_shape_file <- all_shape_files[st_intersects(loc_point,all_shape_files,sparse=F),]
        if(nrow(this_shape_file) > 1){
          warning("The central point for location",location_name,"intersects two shapefiles")
        }
      }
    }
    
    if(method == 'name'){
      ## Get shape file by name
      for(field in names(all_shape_files)[endsWith(names(all_shape_files),suffix = as.character(ISO_A2_level))]){
        tmp_field = paste(field,'tmp',sep='_')
        all_shape_files[[tmp_field]] = all_shape_files[[field]]
        if(mode(field)==mode(character(0))){
          all_shape_files[[tmp_field]] = gsub('[[:punct:]]','',all_shape_files[[tmp_field]])
          all_shape_files[[tmp_field]] = gsub(' ','',all_shape_files[[tmp_field]])
          all_shape_files[[tmp_field]] = iconv(all_shape_files[[tmp_field]],from='UTF-8',to='ASCII//TRANSLIT')
          all_shape_files[[tmp_field]] = toupper(all_shape_files[[tmp_field]])
        }
      }
      this_shape_file = all_shape_files
      if(ISO_A2_level > 0){
        this_shape_file = all_shape_files[
          apply(sapply(names(all_shape_files)[endsWith(names(all_shape_files),suffix='tmp')],function(field){
            eval(parse(text=paste(
              "all_shape_files$",field," == '",ISO_A2[ISO_A2_level],"'",sep=''
            )))
          }),1,function(x){any(x,na.rm=TRUE)}),
        ]
      }
      
      if(nrow(this_shape_file) > 1){
        all_shape_files = this_shape_file
        tmp_names <- names(this_shape_file)[grepl(pattern = "^NAME",names(this_shape_file))]
        tmp_names <- tmp_names[!endsWith(tmp_names,suffix='tmp')] #Remove NAME_ISO_A2_level_tmp
        tmp_names <- tmp_names[!endsWith(tmp_names,suffix='0')] #Remove NAME_0
        tmp_names <- tmp_names[!endsWith(tmp_names,suffix='1')] #Remove NAME_0
  
        index = 1:nrow(this_shape_file) * 0 + 1
        for(name in tmp_names){
          tmp_name = paste(name,"tmp",sep='_')
          all_shape_files[[tmp_name]] = gsub('[[:punct:]]','',all_shape_files[[name]])
          all_shape_files[[tmp_name]] = gsub(' ','',all_shape_files[[tmp_name]])
          all_shape_files[[tmp_name]] = iconv(all_shape_files[[tmp_name]],from='UTF-8',to='ASCII//TRANSLIT')
          all_shape_files[[tmp_name]] = toupper(all_shape_files[[tmp_name]])
          index = index * (all_shape_files[[tmp_name]] %in% location.array)
        }
        index = which(index==1)
        if(length(index) > 1){
          cntry_shapefile <- get_shape_file(paste(location.array[1:3],collapse='_'),taxonomy_dir=taxonomy_dir)
          #' @importFrom sf st_within
          tmp = index[sapply(st_within(this_shape_file,cntry_shapefile),function(x){length(x) > 0})]
          if(length(index) == 0){
            #' @importFrom units set_units 
            #' @importFrom sf st_area
            # @importFrom sf st_difference
            index = index[st_area(st_difference(this_shape_file,cntry_shapefile)) < set_units(1,m^2)]
          }
          if(verbose && (length(index) > 1)){warning("This part not written yet")}
        }
        if(verbose &&(length(index) < 1)){
          warning("This part not written yet")
        }
        this_shape_file = all_shape_files[index,]
      }
      if(nrow(this_shape_file) == 0){
        if(verbose){
          warning("Could not find shapefile for location ",location_name)
        }
      }
  
      this_shape_file = this_shape_file[,!endsWith(names(all_shape_files),'tmp')]
    }
    if(is.null(rc_shape_files)){
      rc_shape_files = this_shape_file
    } else {
      #' @importFrom dplyr intersect
      allowed_names = dplyr::intersect(names(rc_shape_files),names(this_shape_file))
      rc_shape_files = rc_shape_files[,names(rc_shape_files) %in% allowed_names]
      this_shape_file = this_shape_file[,names(this_shape_file) %in% allowed_names]
      rc_shape_files = rbind(rc_shape_files,this_shape_file)
    }

    return(rc_shape_files)
  }
}


#' @name fix_location_name
#' @title fix_location_name
#' @description Attempt to look up a standardized representation of a location.
#' @param location Name of a location as a character
#' @param taxonomy_dir Name of a taxonomy directory (shapefiles will be stored here)
#' @return A character containing the corrected filename or NA in case of failure
fix_location_name = function(location,taxonomy_dir = 'taxonomy-verified',verbose=FALSE){
  if(file.exists(paste(taxonomy_dir,'/',location,'_LOC.csv',sep=''))){
    return(location)
  }
  data('WHO_regions',package='taxdat')

  # location should look like `WHO-REGION_ISO-A1_ISO-A2-L1_...`
  # get each component
  location.array = strsplit(location,split = '_')[[1]]
  if(endsWith(location.array[length(location.array)],'csv')){
    location.array = location.array[-length(location.array)]
  }
  # throw warning if blank
  if(length(location.array) == 0){
    if(verbose){
      warning("Location not provided")
    }
    return(NA)
  }
  #location.array has at least one element

  # starting to check the WHO region
  old.who_region = location.array[1]
  # check to see if WHO region is in long form (example: "AFRICA" instead of "AFR")
  who_region = old.who_region
  if(!(old.who_region %in% who_region_shortener)){
    who_region = who_region_shortener[old.who_region]
  }

  if(verbose){
    if(is.na(old.who_region)){
      message("who region is not a valid region")
    } else if(is.na(who_region)){
      message("who region is not a valid region")
    } else if(old.who_region != who_region){
      message("The who region was changed to the abbreviation.")
    }
  }

  # We're done if its just a WHO region
  if(length(location.array)==1){
    return(who_region)
  }

  #location.array has at least two elements
  #start checking the ISO-A1
  ISO_A1 = location.array[2]
  old.ISO_A1 = ISO_A1

  #Test to see if the ISO_A1 level is part of the ISO_3166_1 standard
  ISO_A1 = toupper(
      gsub('[[:punct:]]',
           '',
           gsub(' ','',
                iconv(ISO_A1,from='UTF-8',
                      to='ASCII//TRANSLIT'))))

  if(any(ISO_A1 == my_ISO_3166_1,na.rm = TRUE)){
    ISO_A1 = ISOcodes::ISO_3166_1[which(ISO_A1 == my_ISO_3166_1,arr.ind = T)[1],'Alpha_3']
  } else if(any(ISO_A1 == my_WHO_regions[,'Entity'],na.rm=TRUE)){
    ISO_A1 = my_WHO_regions[ISO_A1 == my_WHO_regions[,'Entity'],"Country.code"]
  } else {
    if(verbose){
      warning(paste("Could not find the ISO_A1 region",ISO_A1))
    }
    return(NA)
  }
  if(verbose && (old.ISO_A1 != ISO_A1)){
    message("ISO_A1 was not the ISO_3166_1 Alpha 3 code")
  }

  ## We're going to try and pull a Shapefile for more information
  tmp.location = paste(who_region,ISO_A1,sep='_')
  ## the non iso_al info
  ISO_names = location.array[-1]
  #This tests to see if the ISO is valid, but cannot correct it
  # Add code to convert the ISO_A1 region into the country code from various forms including name and ISO2 codes
  err <- tryCatch({
    shapefile = get_shape_file(tmp.location, taxonomy_dir)
    #This will get the standard ISO_A1 region
    ISO_names[1] = shapefile[["NAME_ENGLISH"]]
    0
  }, error = function(e) {
    warning(e$message)
    if(verbose){
      message("The ISO_A1 region", tmp.location, "is not valid.")
    }
    return(1)
  })

  if(err != 0){return(NA)}

  # Now we go back and make sure the who region we have is the same as the one associated with the ISO_A1
  old.who_region = who_region
  who_region = who_region_shortener[
    filter(WHO_regions,Country.code == ISO_A1)$WHO.region
    ]
  if(verbose){
    if(is.na(old.who_region)){
      message("who region is not a valid region")
    } else if(old.who_region != who_region){
      message("The who region was changed to match the country")
    }
  }

  # Return if done
  if(length(location.array) == 2){
    return(paste(who_region,ISO_A1,sep='_'))
  }
  #location.array has at least three elements
  ISO_A2 = location.array[3:length(location.array)]
  old.ISO_A2 = ISO_A2
  ISO_A2_level = length(location.array)-2

  if(!file.exists(paste(taxonomy_dir,"ShapeFiles",sep='/'))){
    stop("The directory",taxonomy_dir,"should exists, and have a subdirectory called ShapeFiles")
  }
  # Now we check each ISO_A2 level to make sure it exists.
  for(level in 1:ISO_A2_level){

    ISO_A2[level] = toupper(
      gsub('[[:punct:]]',
           '',
           gsub(' ','',
                iconv(ISO_A2[level],from='UTF-8',
                      to='ASCII//TRANSLIT'))))

    destination = paste(taxonomy_dir,'/ShapeFiles/',ISO_A1,'_adm',level,'.rds',sep='')
    if(!file.exists(destination)){
      website = paste("http://biogeo.ucdavis.edu/data/gadm2.8/rds/",ISO_A1,"_adm",level,".rds",sep='')
      try({download.file(website,destination,mode='wb')})
    }
    if(!file.exists(destination)){
      if(verbose){
        message(paste("Could not find shape file at ISO_A2 level",level))
      }
      return(NA)
    }
    #' @importFrom sf st_as_sf
    all_shape_files = st_as_sf(readRDS(destination))
    if(!all(suppressWarnings(st_is_valid(all_shape_files)))){
      #' @importFrom lwgeom st_make_valid
      all_shape_files <- suppressWarnings(st_make_valid(all_shape_files))
    }
    relevent_fields = names(all_shape_files)[endsWith(names(all_shape_files),suffix = as.character(level))]
    for(field in relevent_fields){
      tmp_field = paste(field,'tmp',sep='_')
      all_shape_files[[tmp_field]] = all_shape_files[[field]]
      if(mode(all_shape_files[[field]])==mode(character(0))){
        all_shape_files[[tmp_field]] = strsplit(all_shape_files[[tmp_field]],'|',fixed=TRUE)
        all_shape_files[[tmp_field]] = lapply(all_shape_files[[tmp_field]],function(x){gsub(' ','',x)})
        all_shape_files[[tmp_field]] = lapply(all_shape_files[[tmp_field]],function(x){iconv(from='UTF-8',to='ASCII//TRANSLIT',x)})
        all_shape_files[[tmp_field]] = lapply(all_shape_files[[tmp_field]],function(x){gsub('[[:punct:]]','',x)})
        all_shape_files[[tmp_field]] = lapply(all_shape_files[[tmp_field]],function(x){toupper(x)})
        all_shape_files[[tmp_field]][sapply(all_shape_files[[tmp_field]],length) == 0] = ''
      }
    }

    this_shape_file = all_shape_files[apply(
      sapply(
        names(all_shape_files)[endsWith(names(all_shape_files),suffix='tmp')],
        function(field){
          sapply(
            1:length(all_shape_files[[field]]),
            function(index){
              if(mode(all_shape_files[[field]]) == 'list'){
                ISO_A2[level] %in% all_shape_files[[field]][[index]]
              } else {
                ISO_A2[level] == all_shape_files[[field]][index]
              }
            }
          )
        }
      ),1,function(x){any(x,na.rm=TRUE)}),]

    if(nrow(this_shape_file) == 0){
      if((level==1) & (ISO_A2[level] %in% my_ISO_3166_2[,'Code'])){
        old.ISO_A2_L1 = ISO_A2[level]
        ISO_A2[level] = my_ISO_3166_2[(ISO_A2[level] == my_ISO_3166_2[,'Code']),"Name"]
        this_shape_file = all_shape_files[apply(
          sapply(
            names(all_shape_files)[endsWith(names(all_shape_files),suffix='tmp')],
            function(field){
              sapply(
                1:length(all_shape_files[[field]]),
                function(index){
                  if(mode(all_shape_files[[field]]) == 'list'){
                    ISO_A2[level] %in% all_shape_files[[field]][[index]]
                  } else {
                    ISO_A2[level] == all_shape_files[[field]][index]
                  }
                }
              )
            }
          ),1,function(x){any(x,na.rm=TRUE)}),]
        if(nrow(this_shape_file) == 0){
          if(verbose){
            warning(paste("Could not find shapefile corresponding to ISO_A2_L1 region",ISO_A2[level],"."))
          }
          return(NA)
        }

        if(nrow(this_shape_file) > 1){
          if(verbose){
            warning(paste("Found too many shapefiles corresponding to ISO_A2_L1 region",ISO_A2[level],"."))
          }
          # return(NA)
        }

        ISO_A2[level] = old.ISO_A2_L1
        ##For this case, we can look up this ISO_31666_2 code
      } else {
        if(verbose){
          warning(paste("Could not find ISO_A2_L",level," region.",sep=''))
        }
        return(NA)
      }
    }

    ISO_A2[level] = this_shape_file[[paste('NAME',level,sep='_')]]
    ISO_A2[level] = gsub(' ','-',ISO_A2[level])
    ISO_names[level+1] = this_shape_file[[paste('NAME',level,sep='_')]]
  }


  if(verbose && (any(ISO_A2 != old.ISO_A2))){
    message("changed ISO_A2 information")
  }

  tmp.location = paste(who_region,ISO_A1,paste(ISO_A2[1:ISO_A2_level],collapse='_'),sep='_')
  shapefile = as.data.frame(get_shape_file (tmp.location,taxonomy_dir))
  if(nrow(shapefile) == 0){
    if(verbose){
      warning("The location in question could not be found")
    }
    return(NA)
  }
  if(nrow(shapefile) > 1){
    index = which(apply(
      sapply(0:ISO_A2_level,function(x){shapefile[[paste("NAME",x,sep='_')]] %in% ISO_names}),
      1,
      prod
    )==1)
    if(length(index) > 1){
      warning("This part isn't written yet.")
      return(NA)
    }
    if(length(index) < 1){
      warning("This part isn't written yet.")
      return(NA)
    }
    shapefile <- shapefile[index,]
  }
  ISO_A2 = shapefile[,startsWith(names(shapefile),'NAME') & (!endsWith(names(shapefile),'0')) & (!endsWith(names(shapefile),'tmp'))]

  if(verbose && (!all(old.ISO_A2 %in% ISO_A2))){
    message("The ISO regions were changed when looking up the shapefile")
  }

  old.ISO_A2 = ISO_A2
  ## Use this later
  ISO_A2 = old.ISO_A2[1]
  tmp.location = paste(who_region,ISO_A1,paste(ISO_A2,collapse='_'),sep='_')
  shapefile = get_shape_file(tmp.location,taxonomy_dir)
  if(nrow(shapefile) == 0){
    if(verbose){
      warning(paste("Could not find shpaefile corresponding to ISO_A1 region",tmp.location))
    }
    return(NA)
  }
  ISO_A2[1] = gsub('\\.','-',shapefile$HASC_1)
  if(is.na(ISO_A2[1])){
    if(verbose){
      warning("No 2 letter code provided.  Returning NULL")
    }
    return(NA)
  }
  if(verbose && (old.ISO_A2[1] != ISO_A2[1])){
    message("Changed ISO_A2_L1 to abbreviated form.")
  }
  old.ISO_A2[1] = ISO_A2[1]
  ISO_A2 = old.ISO_A2
  location = paste(who_region,ISO_A1,paste(ISO_A2,collapse='_'),sep='_')
  if(nchar(location) < 5){
    # browser()
  }
  location.tmp = location
  location.tmp = gsub('-','dashcharacter',location.tmp)
  location.tmp = gsub('_','underscorecharacter',location.tmp)
  location.tmp = gsub("[[:punct:]]",'',location.tmp)
  location.tmp = gsub(' ','',location.tmp)
  location.tmp = iconv(location.tmp, from="UTF-8",to="ASCII//TRANSLIT")
  location.tmp = gsub('dashcharacter','-',location.tmp)
  location.tmp = gsub('underscorecharacter','_',location.tmp)
  if(verbose && (location != location.tmp)){
    message("Changed location to be ascii and punctuation free")
  }
  return(location.tmp)

  return(location)
}


#' @name exists_shape_file
#' @title exists_shape_file
#' @description Check to see if a custom shapefile exists in the taxonomy already.
#' @param location Name of a location as a character
#' @param taxonomy_dir Name of a taxonomy directory (shapefiles will be stored here)
#' @param date which date to look for a custom shapefile covering.
#' @return boolean TRUE if a custom shapefile exists, FALSE otherwise
#' @importFrom lubridate now
exists_shape_file = function(location,taxonomy_dir = 'taxonomy-verified',verbose=FALSE,date=now()){
  file = paste0(taxonomy_dir,'/Location/',location,'_LOC.csv')
  ## Check if there is a shapefile
  if(!file.exists(file)){
    if(verbose){
      warning(paste("Could not find location file",file))
    }
    return(FALSE)
  }
  loc_data <- read_location_csv(file)
  ## Check if there is a shapefile field
  if(!any(grepl('gis',names(loc_data)))){
    return(FALSE)
  }
  
  if(all(is.na(loc_data[,grepl('gis',names(loc_data))]))){
    return(FALSE)
  }
  
  ## Check for time limitations on the shapefile
  if(!(
    any(grepl('gis_start',names(loc_data))) |
    any(grepl('gis_end',names(loc_data)))
  )){
    return(loc_data['gis_file'])
  }
  starts = as.matrix(loc_data)[1,grepl('gis_start',names(loc_data))]
  ends = as.matrix(loc_data)[1,grepl('gis_end',names(loc_data))]
  files = as.matrix(loc_data)[1,grepl('gis_file',names(loc_data))]
  if(length(starts) != length(ends)){
    stop("This should not happen.  Please make sure that",file,"has the same number of gis_start and gis_end lines.")
  }
  if(length(files) != length(starts)){
    stop("This should not happen.  Please make sure that",file,"has the same number of gis_file and gis_end lines.")
  }
  gis_file = FALSE
  for(i in 1:length(files)){
    if(
      ((starts[paste0('gis_start_',i)] < date) | is.na(starts[paste0('gis_start_',i)])) &
      ((starts[paste0('gis_end_',i)] < date) | is.na(starts[paste0('gis_end_',i)]))
    ){
      if(i==1){
        gis_file = files['gis_file']
      } else {
        gis_file = files[paste0('gis_file_',i)]
      }
    }
  }
  return(gis_file)
}


#' @name read_location_csv
#' @title read_location_csv
#' @description This function reads a csv file containing a location file
##  args:
#' @param filename A string for the relative or absolute path to the
#'   file to be read.  The file format is described in the
#'   documentation
##  vals:
#' @return A \code{data.frame} containing the data from
#'   \code{filename} or a single column missing with the filename
#'   listed
################################################################

location_coltypes = c(
  name = "character",
  cent_lat = "numeric",
  cent_long = "numeric",
  isISO_A1 = "numeric",
  isISO_A2_L1 = "numeric",
  isISO_A2_L2 = "numeric",
  isISO_A2_L3 = "numeric",
  isISO_A2_L4 = "numeric",
  isISO_A2_L5 = "numeric",
  isISO_A2_L6 = "numeric",
  isISO_A2_L7 = "numeric",
  isISO_A2_L8 = "numeric",
  isISO_A2_L9 = "numeric",
  gis_file = 'character',
  gis_file_2 = 'character',
  gis_file_3 = 'character',
  gis_file_4 = 'character',
  gis_file_5 = 'character',
  gis_start_1 = 'Date',
  gis_start_2 = 'Date',
  gis_start_3 = 'Date',
  gis_start_4 = 'Date',
  gis_start_5 = 'Date',
  gis_end_1 = 'Date',
  gis_end_2 = 'Date',
  gis_end_3 = 'Date',
  gis_end_4 = 'Date',
  gis_end_5 = 'Date',
  enclosed_by = "character",
  notes = "character"
)

read_location_csv = function(filename){
  rc <- read_transposed_taxonomy_csv(filename)
  for(column in colnames(rc)){
    if(!is.na(location_coltypes[column])){
      if(location_coltypes[column] == "Date"){
        # rc[[column]] <- as.Date(rc[[column]])
        rc[[column]] <- safe.function(as.Date,x=rc[[column]],default=as.Date(NA))
      } else {
        rc[[column]] <- as(rc[[column]],location_coltypes[column])
      }
    }
  }
  if(nrow(rc) > 1){
    return(data.frame(missing=filename,stringsAsFactors = FALSE))
  }

  old.ncol = ncol(rc)
  if(old.ncol == 1){
    return(rc)
  }
  rc <- rc[,!is.na(colnames(rc))]
  rc <- rc[,!(colnames(rc) == "")]
  if(ncol(rc) != old.ncol){
    warning(paste("Columns were removed from file",filename,"during the reading process.  They are presumed to be empty"))
  }  
  ## rc <- rc %>% mutate(is_public = is_public == 1)
  return(rc)
}


#' @name read_transposed_taxonomy_csv
#' @title read_transposed_taxonomy_csv
#' @description This function reads a csv file which is row major.  It
#'   will also fail without throwing an error, so it can be used on
#'   missing data
##  args:
#' @param filename A string for the relative or absolute path to the
#'   file to be read.  The file should be row major
##  vals:
#' @return A data.frame containing the data from \code{filename} or a 
#'   data.frame with a single column missing=\code{filename})
################################################################

read_transposed_taxonomy_csv = function(filename){
  if(!file.exists(filename)){
    warning(paste(filename,": Does not exist"),immediate.=FALSE)
    return(data.frame(missing=filename,stringsAsFactors = FALSE))
  }
  output.data  = tryCatch(
    as.data.frame(
      t(read.csv(
        filename,
        sep=',',
        header=FALSE,
        stringsAsFactors=FALSE,
        na.strings = "",
        colClasses = 'character',
        quote = "\"",
        check.names = FALSE
      )),
      stringsAsFactors=FALSE
    ),
    warning = function(w){
      if(length(w$message > 0)){
        if(grepl(pattern="ncomplete final line",x=w$message)){
          return(suppressWarnings(
            as.data.frame(
              t(read.csv(
                filename,
                sep=',',
                header=FALSE,
                stringsAsFactors=FALSE,
                na.strings = "",
                colClasses = 'character',
                quote = "\"",
                check.names = FALSE
              )),
              stringsAsFators = FALSE
            )))
        }
      }
      warning(paste(filename,":",w$message),immediate.=TRUE)
      return(w)
    },
    error = function(e){
      warning(paste(filename,":",e),immediate.=TRUE)
      return(e)
    }
  )
  if(length(class(output.data))>0){
    if(class(output.data)[[1]] != "data.frame"){
      warning(paste(
        "Output date is of type",
        class(output.data),
        "instead of data.frame."
      ))
      output.data = data.frame(
        missing=filename,stringsAsFactors = FALSE
      )
    } else{
      colnames(output.data) = as.character(unlist(output.data[1,]))
      output.data = output.data[-1,]
    }
  }
  return(output.data)
}