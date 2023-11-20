# clean out the environment
rm(list = ls())

# load packages
library(tidyverse)       # for general data wrangling & plotting
library(furrr)           # for parallelizing
library(feather)         # for reading in and writing out large files faster

# load helper functions
source(file.path("src", "utils", "project_impacts.R"))

# map where GDP projections at ADM0-level are stored
dir_output_data <- file.path("/net", "cfc", "landclim1", "pwaidelich")
dir_files_adm0 <- file.path(dir_output_data, "fulldistr_stagg", "adm0")

# map where global aggregate results will be stored
dir_files_global <- file.path(dir_output_data, "fulldistr_stagg", "global")

# create the directory if it does not exist yet
if(!dir.exists(dir_files_adm0)) dir.create(dir_files_global)

# load model run overview
df_modelscen <- read_csv(file.path("data", "input", "df_modelscen.csv")) %>%
  mutate(identifier = paste0(model, "_", scenario, "_", ensemble))

# GDP weights for ADM0-to-global aggregation
df_ssp <- readRDS(file.path("data", "df_ssp.rds"))

# initiate parallel processing
plan(multisession(workers = 40))



################################################################################
############### AGGREGATE ADM0 TO GLOBAL  ######################################
################################################################################

future_map(df_modelscen$identifier,
    .f = function(identifier_selected) {
      
      # print the ID to track progress
      print(identifier_selected)
      
      # detect all files in the identifier's subfolder
      files_for_identifier <- file.path(dir_files_adm0, identifier_selected) %>% list.files()
      
      #  we expect 1,000 (MC draws) x 231 (ADM0 territories) per model run - throw a warning & skip otherwise
      if(length(files_for_identifier) != 231*1000) {
        
        warning(paste0(identifier_selected, " features ", length(files_for_identifier), " instead of 231,000 files. Please inspect"))
        return(NULL)
        
      }
      
      # create the target directory if required
      if(!dir.exists(file.path(dir_files_global, identifier_selected))) dir.create(file.path(dir_files_global, identifier_selected))
      
      # extract the MC draws from the filenames
      files_mc_draw <- files_for_identifier %>% str_extract("(?<=_mc)[:digit:]+(?=_n1000)") %>% as.integer()
      
      # abort if we do not get the expected output (1-1000, 231 times)
      if(!all.equal(sort(files_mc_draw), rep(1:1000, 231) %>% sort())) stop("Something went wrong with extracting the MC draw from filenames")
      
      # loop through the MC draws and aggregate them up
      map(1:1000, .f = function(monte_carlo_draw_selected) {
        
        # pick all filenames for the respective MC draw
        files_mc_draw <- files_for_identifier[files_mc_draw == monte_carlo_draw_selected]
        
        # we expect 231 matches (one for each ADM0 region) - abort otherwise
        if(length(files_mc_draw) != 231) stop("files_mc_draw must have length 231")
        
        # extract a string of how to name the global file
        filename_draw <- files_mc_draw %>% str_extract("mc[:digit:]+_n1000\\.feather$") %>% unique()
        
        # abort if we do not have a unique match
        if(length(filename_draw) != 1) stop("filename_draw is non-unique. Please inspect")
        
        # skip if the global file already exists
        if(file.exists(file.path(dir_files_global, identifier_selected, filename_draw))) {
          paste0("Skipping ", identifier_selected, " ", filename_draw, " which already exists") %>% print()
          return(NULL)
        }
        
        # read in the files, paste them together, aggregate and save out
        map_dfr(files_mc_draw, ~ read_feather(file.path(dir_files_adm0, identifier_selected, .x))) %>%
          aggregate_adm0_to_global(adm0_weights_used = df_ssp, is_mc_output = T) %>%
          write_feather(file.path(dir_files_global, identifier_selected, filename_draw))
        
        return(TRUE)
      # end of inner map function
      })
    
      return(TRUE)
}) 
