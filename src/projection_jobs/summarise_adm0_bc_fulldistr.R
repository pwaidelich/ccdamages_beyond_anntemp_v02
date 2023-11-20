## RUNTIME: 2h using 36 workers

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

# read in the list of all ISO3 codes in the projections
sovereign_countries_gid0 <- readRDS(file.path("data", "sovereign_countries_gid0.rds"))

# read in the GCM weights for the distribution incl. large ensembles
df_gcm_weights_le <- readRDS(file.path("data", "df_gcm_weights_le.rds"))

# initiate parallel processing
plan(multisession(workers = 44))



################################################################################
############### IMPORT & SUMMARIZE  ############################################
################################################################################

# map all ADM0-level countries for which full MC distribution results are compiled
gid0_available <- file.path(dir_output_data, "fulldistr_stagg", "adm0_gwl_feather") %>% list.files() %>%
  str_extract("^[A-Z][A-Z][A-Z](?=\\.feather)")

# loop through them and calculate summary statistics
future_map(gid0_available, function(gid0_selected) {
  
  # print status to track progress
  print("Calculating summary stats for:")
  print(gid0_selected)
  
  # skip if the file already exists
  if(file.exists(file.path("data", "df_adm0_bc_fulldistr_le_aggr", paste0(gid0_selected, ".feather")))) {
    print("File already exists. Skipping...")
    return(NULL)
  }
  
  # calculate summary statistics and track computation time
  # NOTE: to speed up computations, we use a lighter version that only calculates the country-level summary
  # stats of interest, namely summarystats_for_distribution_simplified()
  start <- Sys.time()
  
  file.path(dir_output_data, "fulldistr_stagg", "adm0_gwl_feather", paste0(gid0_selected, ".feather")) %>%
    read_feather() %>%
    data.table::as.data.table() %>%
    summarystats_for_distribution_simplified(gcm_weights_used = df_gcm_weights_le,
                                             is_input_datatable = T, warming_levels_used = c(1.5, 2, 3)) %>%
    # save summary stats out as a .feather
    write_feather(file.path("data", "df_adm0_bc_fulldistr_le_aggr", paste0(gid0_selected, ".feather")))
  
  end <- Sys.time()
  
  # print out time required for calculation
  print("Calculation took:")
  print(end - start)
  
  return(TRUE)
})


################################################################################
############### PREPARE ADM0-level SUMMARY STATS  ##############################
################################################################################

# detect all files
filenames_adm0 <- file.path("data", "df_adm0_bc_fulldistr_le_aggr") %>% list.files() %>%
  str_subset("\\.feather$")

# check their length
if(length(filenames_adm0) != length(sovereign_countries_gid0)) stop("No. of ADM0-level files not matching")

# write them out as one feather file
df_adm0_bc_fulldistr_le_aggr <- file.path("data", "df_adm0_bc_fulldistr_le_aggr", filenames_adm0) %>%
  map_dfr(read_feather)

# reformat and convert agreement_negative (= share of distribution predicting negative impacts)
# into agreement_meansign (= share of distribution agreeing with the sign of the mean impact)
df_adm0_bc_fulldistr_le_aggr <- df_adm0_bc_fulldistr_le_aggr %>% setNames(names(.) %>% str_replace("agreement_negative", "agreementnegative")) %>%
  pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "imp") %>%
  mutate(summarystat = str_extract(variable, "(?<=_)[:alnum:]+$"),
         variable = str_remove(variable, "_[:alnum:]+$")) %>%
  pivot_wider(names_from = "summarystat", values_from = "imp") %>%
  mutate(agreement_meansign = if_else(mean < 0, agreementnegative, 1 - agreementnegative)) %>%
  select(-agreementnegative)

# export summary stats for all countries as one .feather file
write_feather(df_adm0_bc_fulldistr_le_aggr, file.path("data", "df_adm0_bc_fulldistr_le_aggr.feather"))

# clean up the environment
rm(df_adm0_bc_fulldistr_le_aggr, filenames_adm0)

