# clean out the environment
rm(list = ls())

# load packages
library(dtplyr)
library(data.table)
library(tidyverse)       # for general data wrangling & plotting
library(furrr)           # for parallelizing
library(feather)         # for reading in and writing out large files faster

# map where GDP projection files (ADM0-level) will be stored
dir_output_data <- file.path("/net", "cfc", "landclim1", "pwaidelich")

# load helper functions
source(file.path("src", "utils", "analyse_impacts.R"))

# load global warming level data for model-scenario-specific baseline periods based on global warming level (defaulting to 0.84 in src/utils/project_impacts() )
df_gwl <- read_csv(file.path("data", "input", "df_gwl.csv"))

# read in the list of all ISO3 codes in the projections
sovereign_countries_gid0 <- readRDS(file.path("data", "sovereign_countries_gid0.rds"))

# map subfolders (= identifiers of model runs used)
subfolders <- file.path(dir_output_data, "fulldistr_stagg", "adm0") %>% list.files()

# initiate parallel processing
plan(multisession(workers = 40))


################################################################################
############### IMPORT & SUMMARIZE  ############################################
################################################################################

# read in files already calculated
gid0_done <- file.path(dir_output_data, "fulldistr_stagg", "adm0_gwl_feather") %>% list.files() %>%
  str_remove("\\.feather$")

# map missing files
gid0_missing <- sovereign_countries_gid0[!sovereign_countries_gid0 %in% gid0_done]

# loop through the rest
map(gid0_missing, .f = function(gid_0_selected) {
   
  # track progress
  print("Compiling results incl. Monte Carlo draws for dose-response function for country:")
  print(gid_0_selected)
  
  # skip if file exists
  if(file.exists(file.path(dir_output_data, "fulldistr_stagg", "adm0_gwl_feather", paste0(gid_0_selected, ".feather")))) {
    print("Skipping as file already exists")
    return(NULL)
  }
  
  # track calculation time for the current ADM0 country
  start <- Sys.time()
  
  tryCatch({
    # step 1: loop through all subfolders (= identifiers of model runs used), grab files belonging to gid_0_selected and compile
    # to speed up file loading, we use parallelization here via furrr
    future_map_dfr(subfolders, ~ {
                     tryCatch({
                       
                       # create the name of all files of interest in the subfolder
                       file_paths <- file.path(dir_output_data, "fulldistr_stagg", "adm0", .x, 
                                               paste0(gid_0_selected, "_mc", 1:1000, "_n1000.feather"))
                       
                       # read in all files and bind them together
                       purrr::map_dfr(file_paths, read_feather)
                       
                     }, error = function(e) {
                       
                       message(paste0("Error processing ", gid_0_selected, e$message))
                       return(tibble())  # return an empty tibble in case of error
                     })
                   }) %>%
    
    # convert to data.table to speed up subsequent calculations
    data.table::as.data.table() %>% dtplyr::lazy_dt() %>%
    
    # merge in GWLs and subset to years in GWL windows
    left_join(data.table::as.data.table(df_gwl), by = c("model", "scenario", "ensemble", "year")) %>%
    filter(warming_level >= 1) %>%
    
    # rename CESM2-LE to CESM2 to ensure that we treat the single run for different RCP-SSP as the same model
    mutate(model = if_else(model == "CESM2-LE", "CESM2", model)) %>%
    
    # collect data.table results and write out
    as_tibble() %>%
    write_feather(file.path(dir_output_data, "fulldistr_stagg", "adm0_gwl_feather", paste0(gid_0_selected, ".feather")))

  }, error = function(e) message(paste0("Error processing", gid_0_selected, e$message)))
  
  # print out calculation time for the current ADM0 country
  end <- Sys.time()
  print("Calculating took:")
  print(end - start)
    
  return(TRUE)
})