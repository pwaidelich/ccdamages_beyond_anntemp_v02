rm(list = ls())

library(stagg)
library(feather)
library(furrr)
library(tidyverse)
library(ncdf4)
library(data.table)
library(dtplyr)
library(ncdf4.helpers) # for extracting time dimension values from a netCDF file with varying calendar formats

# map input and output directory
dir_climate_data_adm1 <- file.path("/net", "cfc", "landclim1", "pwaidelich", "cmip6_ng_adm1_stagg")
dir_climate_data_g025 <- '/net/argon/landclim/fbatibeniz/clim_indices/data/cmip6-ng'

# load the first year of interest based on +0.84C GWL
first_year_GWL084 <- readRDS(file.path("data", "first_year_GWL084.rds"))

# source helper functions for kelvin-to-celsius conversion
source(file.path("src", "utils", "extract_climate_indices.R"))

# load tibble storing all model runs used
df_modelscen <- read_csv(file.path("data", "input", "df_modelscen.csv")) %>%
  # create a model-scenario-ensemble identifier used to name subfolders
  mutate(identifier = paste0(model, "_", scenario, "_", ensemble))

# load the grid cell area weights and adjust the grid from 0-360 to -180-to-180 longitude
df_adm1_weights_stagg <- file.path("data", "df_adm1_weights_stagg.feather") %>% read_feather() %>%
  dplyr::select(lon = "x", lat = "y", GID_1 = "poly_id", w_area) %>%
  mutate(lon = if_else(lon >= 180, lon - 360, lon))

# extract unique latitudes that feature in *any* ADM1 polygon (we discard remaining ones below for efficiency)
lat_with_weights <- unique(df_adm1_weights_stagg$lat) %>% sort()

# write a function for area-weighted aggregation
aggregate_g025_to_adm1 <- function(filepath = NULL,
                                   dir_output = NULL,
                                   varname_kotz_selected = NULL,
                                   identifier_selected = NULL,
                                   cutoff_year = first_year_GWL084,
                                   mask_lat = lat_with_weights) {
  
  # determine the Kotz et al (2022)-style variable name based on RegEx in the filepath
  if(is.null(varname_kotz_selected)) {
    
    varname_kotz <- case_when(str_detect(filepath, "ann_totd_") ~ "Pt",
                              str_detect(filepath, "ann_mean_") ~ "Tmean",
                              # NOTE: not all netCDF files are consistetly named, so we allow for different spellings below
                              str_detect(filepath, "dtd_var_") ~ "Tstd",
                              str_detect(filepath, "dtd_temp_") ~ "Tstd",
                              str_detect(filepath, "btstrp_pr999_") ~ "vwet_days1_am_99p9",
                              str_detect(filepath, "mon_totd_") ~ "mon_totd",
                              str_detect(filepath, "mon_sum_") ~ "mon_totd",
                              str_detect(filepath, "wet_days_") ~ "wet_days_1",
                              str_detect(filepath, "wetd_") ~ "wet_days_1",
                              str_detect(filepath, "[Tt]x5day") ~ "tx5d",
                              TRUE ~ NA_character_
    )
  } else {
    
    # if varname is user-specified, we use that one instead
    varname_kotz <- varname_kotz_selected
  }
  
  # if varname_kotz is still NA at this stage, throw an error
  if(is.na(varname_kotz)) stop("varname_kotz could not be identified in filepath")
  
  # if there is no user-specified model-scenario-ensemble identifier, extract it from filepath using RegEx
  if(is.null(identifier_selected)) {
    
    identifier <- stringr::str_extract(filepath, "(?<=tas_|pr_|tasmax_).+(?=_g025\\.nc$)")
    
    # some files have "_native.nc" instead of "_g025.nc" as an ending, so we allow for this, too
    if(is.na(identifier)) identifier <- stringr::str_extract(filepath, "(?<=tas_|pr_|tasmax_).+(?=_native\\.nc$)")
    
    # if there is still no match, throw an error
    if(is.na(identifier)) stop("Identifier string could not be extracted from filepath")
    
    # if the identifier is *not* included in df_modelscen (= overview of model runs considered), we skip this model-scenario-ensemble pairing
    # NOTE: this serves to prevent aggregating files for CMIP6 models, for which we only have some of the indicators available
    if(!identifier %in% df_modelscen$identifier) {
      
      print(paste0(identifier, " not included in df_modelscen. Skipping..."))
      
      return(NULL)
    }
    
  } else {
    
    # if the user provides an identifier manually (e.g., for ERA5 files), we use that one instead
    identifier <- identifier_selected
  }
  
  # check if a directory named after 'identifier' exists and if not create
  if(!dir.exists(file.path(dir_output, identifier))) {
    
    dir.create(file.path(dir_output, identifier))
  
  } else {
    
    # if the file to be created already exists, we abort to avoid redundant calculations
    if(file.exists(file.path(dir_output, identifier, paste0(varname_kotz, ".feather")))) {
      
      print("File already exists. Skipping...")
      return(NULL) 
    } 
  }
  
  # extract model, scenario & ensemble from the identifier
  # NOTE: this will produce NAs for ERA5 files, which does not cause issues
  model <- str_extract(identifier, "^[^_]+")
  scenario <- str_extract(identifier, "ssp[:digit:][:digit:][:digit:]")
  ensemble <- str_extract(identifier, "(?<=ssp[:digit:][:digit:][:digit:]_)r[:digit:]+i[:digit:]+p[:digit:]+f[:digit:]+")
  
  # open the NetCDF file and extract variable names
  ncfile <- ncdf4::nc_open(filepath)
  varname <- names(ncfile$var)
  
  # get rid of any potential variables that are not tas, pr or tasmax
  varname <- varname[varname %in% c("tas", "pr", "tasmax")]
  
  # check if we have a match
  if(!varname %in% c("tas", "pr", "tasmax")) stop("varname other than tas, pr or tasmax detected. Please inspect the netCDF file")
  
  # get dimensions info for the variable
  dim_info <- ncfile$var[[which(names(ncfile$var) == varname)]]$dim
  
  # check dimension order and lengths
  if(all.equal(c(dim_info[[1]]$name, dim_info[[2]]$name, dim_info[[3]]$name), c("lon", "lat", "time")) != T) stop("Order of dimensions unexpected. Aborting...")
  if(dim_info[[1]]$len != 1440) stop("lon dimension length unexpected")
  if(dim_info[[2]]$len != 721) stop("lon dimension length unexpected")
  
  # check time dimension length
  # NOTE: we have either 1850-2100 (for some files in CESM2-LE 2101) or 1979-2019 (for ERA5)
  # since monthly precip totals have 12x more timesteps, we also allow for each value multiplied by 12
  if(!dim_info[[3]]$len %in% c(251, 41, 251*12, 41*12, 252, 252*12)) {
    
    warning(paste0("time dimension length unexpected. Dim length is ", dim_info[[3]]$len, ". Returning NULL"))
    return(NULL)
  } 
  
  # get the longitude and latitude values
  lon <- dim_info[[1]]$vals
  lat <- dim_info[[2]]$vals
  
  # extract year range
  # NOTE: for 360-day calendars, this throws a warning but result is sound
  time_range <- ncdf4.helpers::nc.get.time.series(ncfile)
  time_subset <- which(lubridate::year(time_range) >= cutoff_year)

  # identify all latitudes with a weight for at least one ADM1 region
  # if our mask features any values that are not actually in the netCDF file, throw an error
  if(mean(mask_lat %in% lat) != 1) stop("Not all values in mask_lat are in the netCDF file. Aborting...")
  lat_subset <- which(lat >= min(mask_lat) & lat <= max(mask_lat))

  # ncdf4::ncvar_get() can specific a start value and a count (see documentation)
  # we set start to the lowest latitude in our mask and the user-specific cutoff-year
  start_nc <- c(1, min(lat_subset), min(time_subset))
  # for count, we use all values for latitude and time but only the length of our latitude mask for latitude
  count_nc <- c(-1, length(lat_subset), -1)
  
  # read the variable of interest using the start and count arguments defined above
  values <- ncdf4::ncvar_get(ncfile, varname, start = start_nc, count = count_nc)

  # close the NetCDF file
  ncdf4::nc_close(ncfile)
  
  # # Filter out times before the cutoff year (e.g., 1964)
  time_range <- time_range[time_subset]
  lat <- lat[lat_subset]

  # for mon_totd we extract year-months, otherwise we extract only years
  if(varname_kotz == "mon_totd") {
    
    year_range <- paste(lubridate::year(time_range), str_pad(lubridate::month(time_range), width = 2, pad = "0"), sep = "-")
    
  } else {
    
    year_range <- lubridate::year(time_range)
    
    # if the extracted year range is neither 1979-2019 (= ERA5) nor cutoff to 2100 (or 2101), throw an error
    if(all.equal(year_range, (cutoff_year):2100) != TRUE & all.equal(year_range, 1979:2019) != TRUE & all.equal(year_range, (cutoff_year):2101) != TRUE)
      stop("Extracted year_range has unexpected values. Please inspect")
  }
  
  ### calculate area-weighted average by GID_1 polygon & year
  # NOTE: we do everything in one step and save out directly to reduce memory loads for parallelization
  
  # create an empty data set with the dimensions of the netCDF using CJ() from data.table()
  data.table::CJ(lon = lon, lat = lat, year = year_range, sorted = F) %>%
    
    # initiate lazy_dt
    dtplyr::lazy_dt() %>%
    
    # re-arrange the output of CJ to replicate the dimension order of the values array extracted from the netCDF file
    arrange(year, lat, lon) %>%
    
    # fill with the values from the netCDF
    dplyr::mutate(value = as.vector(values)) %>%
    
    # merge in grid cell area weights created in a previous script
    dplyr::left_join(data.table::as.data.table(df_adm1_weights_stagg), by = c("lon", "lat")) %>%
    
    # discard grid cells outside of any polygons
    dplyr::filter(!is.na(w_area)) %>%
    
    # calculate weighted average for each ADM1 region and year (NOTE: weights sum to one)
    dplyr::group_by(GID_1, year) %>%
    dplyr::summarise(value = sum(value*w_area)) %>% dplyr::ungroup() %>%
    
    # convert temperature variables to Celsius
    {if(varname_kotz %in% c("Tmean", "tx5d")) dplyr::mutate(., value = kelvin_to_celsius(value)) else . } %>%
    
    # rename the variable to the value of varname_kotz
    dplyr::rename(!!varname_kotz := "value") %>%
    
    # collect the tibble
    dplyr::as_tibble() %>%
    
    # write out as a feather into the subfolder named after identifier and with the varname as filename
    feather::write_feather(file.path(dir_output, identifier, paste0(varname_kotz, ".feather")))
  
  # return TRUE
  return(TRUE)
}


################################################################################
######################## AGGREGATE ERA5 1979-2019 ##############################
################################################################################

# map all relevant files for ERA5
print("Aggregating ERA5 files...")

era5_files <- file.path(dir_climate_data_g025, "era5", "indices") %>%
  list.files(full.names = T) %>%
  # subset to +0.84GWL files
  str_subset("084\\.nc") %>%
  # discard any helper files not used here
  str_subset("stan_mon_dev", negate = T) %>%
  str_subset("threshold", negate = T) %>%
  str_subset("last_era5_tx5day", negate = T)

# apply the function to them
map(era5_files, ~ aggregate_g025_to_adm1(filepath = .x,
                                         dir_output = dir_climate_data_adm1,
                                         varname_kotz_selected = NULL,
                                         identifier_selected = "era5_084"))



################################################################################
######################## AGGREGATE BIAS-CORRECTED CMIP6 OUTPUTS ################
################################################################################

# map all relevant bias-corrected files for CMIP6
all_files <- c( # CMIP6 next-gen single model runs
               file.path(dir_climate_data_g025, "tas", "bc_ann_mean_regrid") %>% list.files(full.names = T),
               file.path(dir_climate_data_g025, "tas", "bc_dtd_var_regrid") %>% list.files(full.names = T),
               file.path(dir_climate_data_g025, "pr", "bc_ann_totd_regrid") %>% list.files(full.names = T),
               file.path(dir_climate_data_g025, "pr", "bc_wet_days_regrid") %>% list.files(full.names = T),
               file.path(dir_climate_data_g025, "pr", "bc_mon_totd_regrid") %>% list.files(full.names = T),
               file.path(dir_climate_data_g025, "pr", "bc_btstrp_pr999_regrid") %>% list.files(full.names = T),
               # MPI large ensemble
               file.path(dir_climate_data_g025, "tas", "bc_ann_mean_mpi_regrid") %>% list.files(full.names = T),
              file.path(dir_climate_data_g025, "tas", "bc_dtd_var_mpi_regrid") %>% list.files(full.names = T),
              file.path(dir_climate_data_g025, "pr", "bc_ann_totd_mpi_regrid") %>% list.files(full.names = T),
              file.path(dir_climate_data_g025, "pr", "bc_wet_days_mpi_regrid") %>% list.files(full.names = T),
              file.path(dir_climate_data_g025, "pr", "bc_mon_totd_mpi_regrid") %>% list.files(full.names = T),
              file.path(dir_climate_data_g025, "pr", "bc_btstrp_pr999_mpi_regrid") %>% list.files(full.names = T),
              # CESM2-LE large ensemble
              file.path(dir_climate_data_g025, "tas", "bc_ann_mean_ensemble_regrid") %>% list.files(full.names = T),
              file.path(dir_climate_data_g025, "tas", "bc_dtd_var_ensemble_regrid") %>% list.files(full.names = T),
               file.path(dir_climate_data_g025, "pr", "bc_ann_totd_ensemble_regrid") %>% list.files(full.names = T),
              file.path(dir_climate_data_g025, "pr", "bc_wet_days_ensemble_regrid") %>% list.files(full.names = T),
              file.path(dir_climate_data_g025, "pr", "bc_mon_totd_ensemble_regrid") %>% list.files(full.names = T),
              file.path(dir_climate_data_g025, "pr", "bc_btstrp_pr999_ensemble_regrid") %>% list.files(full.names = T)
               )

# print a status update
print("Aggregating CMIP6 files w/o large ensembles...")

# loop through rows of df_modelscen
# NOTE: monthly precip totals files take considerably longer and are memory-intensive, so scope for parallelization here is limited
map(1:nrow(df_modelscen),
           .f = function(row_number) {
             
             # print status update
             paste0("Aggregating ", df_modelscen$identifier[row_number], "...") %>% print()

             # map all files for the model-scenario-ensemble combination in the respective row of df_modelscen
             files_g025 <- all_files %>% str_subset(paste0(df_modelscen$model[row_number], "_", df_modelscen$scenario[row_number], "_", df_modelscen$ensemble[row_number]))

             # for first realization of MPI ensemble, we get two results in different folders, so we discard one of them (to avoid finding 12 files instead of 6)
             if(df_modelscen$identifier[row_number] == "MPI-ESM1-2-LR_ssp370_r1i1p1f1") {
               files_g025 <- files_g025 %>% str_subset("mpi_regrid")
             }

             # abort and throw a warning if we do not get exactly 6 matches
             if(length(files_g025) != 6) {
               warning(paste0("Getting more or less than six matches for ", df_modelscen$identifier[row_number], ". Skipping..."))
               return(NULL)
             }

             # apply the aggregate function
             map(files_g025, ~ aggregate_g025_to_adm1(filepath = .x,
                                                      dir_output = dir_climate_data_adm1,
                                                      varname_kotz_selected = NULL,
                                                      identifier_selected = NULL))

            # free up memory via garbage collector
            gc()
            
            # return TRUE
            return(TRUE)
})



################################################################################
######################## AGGREGATE INDIVIDUAL FILES IN CASE OF ERRORS ##########
################################################################################

# NOTE: uncomment this and run with user-specific values for identifier
# if individual files were not aggregated for some reason as part of the map() loop above

### UNCOMMENT CODE BELOW
# identifier <- "MPI-ESM1-2-LR_ssp370_r2i1p1f1"
# 
# # catch some individual manually if required (e.g., because feather files are corrupted)
# aggregate_g025_to_adm1(filepath = file.path(dir_climate_data_g025, "tas", "bc_ann_mean_mpi_regrid") %>% list.files(full.names = T) %>% str_subset(identifier),
#                        dir_output = dir_climate_data_adm1,
#                        varname_kotz_selected = NULL,
#                        identifier_selected = NULL)

### END OF CODE TO UNCOMMENT


################################################################################
######################## AGGREGATE RAW CMIP6 MONTHLY PRECIP TOTALS #############
################################################################################

# parallelize across two workers to speed this up
plan(multisession(workers = 2))

# map all relevant files
all_files_Wn_raw <- c(file.path(dir_climate_data_g025, "pr", "mon_totd_ensemble_regrid") %>% list.files(full.names = T),
                      file.path(dir_climate_data_g025, "pr", "mon_totd_mpi_regrid") %>% list.files(full.names = T),
                      file.path(dir_climate_data_g025, "pr", "mon_totd_regrid") %>% list.files(full.names = T)
                      )

# loop through them
future_map(all_files_Wn_raw,
    .f = function(filepath) {
      
      # print filepath to track progress
      print(filepath)
      
      # aggregate
      aggregate_g025_to_adm1(filepath = filepath,
                             # NOTE: we save raw CMIP6 files into a different folder
                             dir_output = file.path(dir_climate_data_adm1 %>% str_replace("stagg$", "stagg_raw")),
                             varname_kotz_selected = NULL,
                             identifier_selected = NULL)
    })



################################################################################
######################## AGGREGATE BIAS-CORRECTED Tx5d #########################
################################################################################

# apply the same loop to all bias-corrected Tx5d netCDF files (variable name in netCDF: tasmax)
map(file.path(dir_climate_data_g025, "tasmax", "bc_Tx5day_regrid") %>% list.files(full.names = T),
    .f = function(filepath) {

      print(filepath)

      aggregate_g025_to_adm1(filepath = filepath,
                             # NOTE: we save Tx5d files into the regular directory, along with our 6 main indicators
                             dir_output = dir_climate_data_adm1,
                             varname_kotz_selected = NULL,
                             identifier_selected = NULL)
    })



################################################################################
######################## AGGREGATE MPI RUN W/ ALTERNATIVE BIAS CORRECTION ######
################################################################################

## apply to example model run with bias correction at daily precip/temp level
map(c(list.files('/net/argon/landclim/fbatibeniz/clim_indices/data/cmip6-ng/tas/indices_bias_cor_var',
               full.names = T) %>% str_subset("MPI"),
      list.files('/net/argon/landclim/fbatibeniz/clim_indices/data/cmip6-ng/pr/indices_bias_cor_var',
                 full.names = T) %>% str_subset("MPI")),
    
    .f = function(filepath) {

      print(filepath)

      aggregate_g025_to_adm1(filepath = filepath,
                             # NOTE: we save alternative bias correction values into a different folder
                             dir_output = dir_climate_data_adm1 %>% str_replace("_stagg", "_stagg_altbiascor"),
                             varname_kotz_selected = NULL,
                             identifier_selected = "MPI-ESM1-2-LR_ssp370_r1i1p1f1")
    })
