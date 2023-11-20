library(tidyverse)  # general data wrangling and plotting
library(ncdf4)      # opening and inspecting netCDF files
library(metR)       # converting netCDF fast and conveniently into a data.frame
library(furrr)      # for parallelizing the map() function
library(ncdf4.helpers) # for extracting time dimension values from a netCDF file with varying calendar formats

# lubridate package is used w/o loading it to avoid namespace issues - throw an error if user does not have it installed
if(!"lubridate" %in% rownames(installed.packages())) {
  stop("Using the extract_adm1_var() function requires that the lubridate package is installed. Please install it and try again")
}

# NOTE: the functions in this file match filenames based on model and scenario regex patterns. This implicitly assumes that
# the directory with ADM1-level climate indices contains only one ensemble ID per model-scenario combination. If this is not
# the case, the functions might produce wrong results although some basic error-checking regarding multiple matches is carried out

# write a function to convert temperature in Kelvin to Celsius
kelvin_to_celsius <- function(x) x - 273.15

# write a easier wrapper that simply collects stagg output
collect_climdata_feathers <- function(model = NULL, scenario = NULL, ensemble = NULL,
                                      dir_files = NULL, identifier_selected = NULL,
                                      include_tx5d = FALSE,
                                      simulate_Wn = FALSE) {
  
  # create identifier
  if(!is.null(identifier_selected)) {
    identifier <- identifier_selected
  } else {
    identifier <- paste0(model, "_", scenario, "_", ensemble)
  }
  
  # read in the filepaths
  files <- file.path(dir_files, identifier) %>% list.files(full.names = T) %>% str_subset("mon_totd", negate = T)
  
  # if we simulate monhtly precip deviation, we discard the respective feather file
  if(simulate_Wn) files <- files %>% str_subset("Wn\\.feather", negate = T)
  
  # unless we include tx5d, we discard this file as well
  if(!include_tx5d) files <- files %>% str_subset("[Tt]x5d", negate = T)
  
  # check for the number of remaining files (should be 6, with tx5d 7, one less if we simulate Wn)
  if(length(files) != 6 - simulate_Wn  & !include_tx5d) stop(paste0("No. of files unexpected for ", identifier))
  if(length(files) != 7 - simulate_Wn & include_tx5d) stop(paste0("No. of files unexpected for ", identifier))
  
  # read in all files, merge together, add identifier info
  df_out <- files %>% map(read_feather) %>% reduce(left_join, by = c("year", "GID_1")) %>%
    mutate(model = model, scenario = scenario, ensemble = ensemble) %>%
    # filter out years over 2100 (since the CESM2-LE has 2101 for some variables)
    filter(year <= 2100)
  
  # simulate monthly precip deviation if specified
  if(simulate_Wn) df_out <- df_out %>% mutate(Wn = 0)
  
  #return
  return(df_out %>% 
           select(model, scenario, ensemble, GID_1, year,
                  Tmean, Tstd, Pt, Wn, wet_days_1, vwet_days1_am_99p9, everything()))
  
}


### PREVIOUS FUNCTIONS - NO LONGER USED AS OF REVISION #1

# extract a selected climate index for a given model and scenario (based on James' code)
extract_adm1_var <- function(varname = NULL, dir = NULL,
                             model = NULL, vars_object = NULL, scenario = NULL,
                             ensemble = NULL) {
  if(is.null(varname) | is.null(dir) | is.null(model) | is.null(vars_object) | is.null(scenario)) stop("Some function arguments have not been specified. Please revise")
  
  # extract the corresponding names in netcdf4 variables and the Kotz et al regression model from the vars object
  var_short <- vars_object$short[vars_object$names == varname]
  name_kotz <- vars_object$names_kotz[vars_object$names == varname]
  
  # open the netcdf file and extract the relevant information
  filename <- file.path(dir, varname) %>% list.files(full.names = TRUE) %>% str_subset("\\.nc$") %>%
    str_subset(paste0("_", model, "_")) %>% str_subset(paste0("_", scenario, "_"))
  
  # subset further if ensemble is specified
  if(!is.null(ensemble)) filename <- filename %>% str_subset(paste0("_", ensemble, "_"))
  
  # ensure the file exists
  if(!file.exists(filename)) stop(paste0("The file ", filename, " does not exist. Please inspect"))
  
  # check if there is one unique file that was matched and throw an error if not
  if(length(filename) == 0) stop(paste("No file was detected for", varname, "from", model, scenario, "in", dir))
  if(length(filename) > 1) stop(paste("Multiple files were detected for", varname, "from", model, scenario, "in", dir))
  
  # open the netCDF file and extract the variable of interest
  nc <- ncdf4::nc_open(filename) 
  values <- ncdf4::ncvar_get(nc, var_short)
  
  # for monthly precip deviation, the time dimension is already in years. For raw CMIP6 files, we extract years by converting from 'days since' to year
  if(name_kotz == "Wn") {
    
    year_range <- ncvar_get(nc, "year")
    
  } else {
    
    year_range <- lubridate::year(nc.get.time.series(nc))
  
  }
  
  # extract the ADM1 region identifier
  GID_1 <- ncdf4::ncvar_get(nc, 'GID_1')
  
  # combine all relevant info into a long-format tibble and return it
  out <- values %>% as_tibble() %>%
    setNames(GID_1) %>%
    mutate(year = year_range) %>%
    gather("GID_1", "varname", -c(year)) %>%
    setNames(names(.) %>% str_replace("varname", name_kotz)) %>%
    mutate(model = model, scenario = scenario)
  
  if(!is.null(ensemble)) out <- out %>% mutate(ensemble = ensemble)
  
  # convert annual average temperature from Kelvin to Celsius
  if(name_kotz == "Tmean") out <- out %>% mutate(Tmean = kelvin_to_celsius(Tmean))
  if(name_kotz == "tx5d") out <- out %>% mutate(tx5d = kelvin_to_celsius(tx5d))
  
  return(out)
}


# write a wrapper that does this for all climate indices of a given model
extract_adm1_allvars <- function(model = NULL, scenario = NULL, dir = NULL,
                                 vars_object = NULL, ensemble = NULL) {
  if(is.null(model) | is.null(scenario) | is.null(vars_object) | is.null(dir)) stop("Please specify the arguments model, scenario, vars_object and dir")
  
  # NOTE: we do NOT parallelize at this stage (across different climate indices) - parallelization occurs at the highest level possible
  out <- map(vars_object$names, ~ extract_adm1_var(varname = .x, dir = dir, vars_object = vars_object,
                                          model = model, scenario = scenario, ensemble = ensemble)) %>%
    # combine the list produced by map() into a single data frame via a left_join
    reduce(left_join, by = c("year", "GID_1", "model", "scenario")) %>%
    {if(!is.null(ensemble)) select(., -starts_with("ensemble")) else .}

  if(!is.null(ensemble)) out <- out %>% mutate(ensemble = ensemble)
    
  return(out %>%
         # reorder the columns
         select(model, scenario, GID_1, year, everything()))
}
