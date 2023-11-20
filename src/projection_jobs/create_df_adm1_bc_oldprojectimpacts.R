## RUNTIME: 6min at 40 workers

# clean out the environment
rm(list = ls())

# load packages
library(tidyverse)       # for general data wrangling & plotting
library(furrr)           # for parallelizing
library(feather)         # for reading in and writing out large files faster

# load helper functions
source(file.path("src", "utils", "extract_climate_indices.R"))
source(file.path("src", "utils", "project_impacts.R"))
source(file.path("src", "utils", "analyse_impacts.R"))

# load global warming level data for model-scenario-specific baseline periods based on global warming level (defaulting to 0.84 in src/utils/project_impacts() )
df_gwl <- read_csv(file.path("data", "input", "df_gwl.csv"))

# GDP weights for ADM1-to-ADM0 aggregation
gdp_weights <- read_csv(file.path("data", "adm1-weights.csv"))  %>%
  group_by(GID_0) %>% mutate(adm1_weight = GDP/sum(GDP)) %>% ungroup() %>%
  select(GID_0, GID_1, GDP, adm1_weight)

# # GDP weights for ADM0-to-global aggregation
df_ssp <- readRDS(file.path("data", "df_ssp.rds"))

# draw coefficients
# NOTE: get_coefficients() uses a default seed for Monte Carlo draws from multi-variate Gaussian, so we do not need to fix the seed here)
coefs_point <- get_coefficients()
coefs_point_T_Pt_Pt2zero <- get_coefficients(model_used = "T_Pt_Pt2zero")

# map the directory where climate projection files are stored
dir_climate_data_adm1 <- file.path("/net", "cfc", "landclim1", "pwaidelich", "cmip6_ng_adm1_stagg")

# map the output directory where GDP projection files (ADM0-level) will be stored
dir_output_data <- file.path("/net", "cfc", "landclim1", "pwaidelich")

# load the tibble objects storing model-scenarios for which we have all 6 climate indicators available
df_modelscen <- read_csv(file.path("data", "input", "df_modelscen.csv"))

# read in the cap for monthly precip deviation from src/00 Create Projections.R
ceiling_Wn_bc <- readRDS(file.path("data", "ceiling_Wn_bc.rds"))

# initiate parallel processing
plan(multisession(workers = 40, gc = TRUE))



################################################################################
############### CREATE PROJECTIONS #############################################
################################################################################

# ensure target directories exist
dir_adm1 <- file.path(dir_output_data, "pointestimates_stagg_oldprojectimpacts", "adm1")
if(!dir.exists(dir_adm1)) dir.create(dir_adm1) 

dir_target_adm0 <- file.path(dir_output_data, "pointestimates_stagg_oldprojectimpacts", "adm0")
if(!dir.exists(dir_target_adm0)) dir.create(dir_target_adm0)  

# loop through all models
future_pmap(.l = df_modelscen,
            .f = function(model, scenario, ensemble) {
              
              # create identifier character
              identifier <- paste0(model, "_", scenario, "_", ensemble)
              
              # print out progress - commented out since this only works for map(), not future_map()
              # paste0("Calculating point estimate results for ", identifier) %>% print()
              
              # calculate ADM1-level impacts using point estimates
              df_adm1 <- collect_climdata_feathers(model = model,
                                                   scenario = scenario,
                                                   ensemble = ensemble,
                                                   dir_files = dir_climate_data_adm1) %>%
                # apply the ceiling for bias-corrected monthly precip deviation
                mutate(Wn = if_else(Wn > ceiling_Wn_bc, ceiling_Wn_bc, Wn)) %>%
                # calculate GDP impacts
                project_impacts_initial_submission(coefs = coefs_point[, 1],
                                return_helper_columns = F,
                                df_gwl_baseline = df_gwl,
                                use_gwl_baseline = T
                )
              
              # determine the target directory for ADM1-level results
              dir_target_adm1 <- file.path(dir_adm1, identifier)
              if(!dir.exists(dir_target_adm1)) dir.create(dir_target_adm1)
              
              # write out in ADM0-level chunks (for easy access)
              df_adm1 %>%
                mutate(GID_0 = str_extract(GID_1, "^[A-Z][A-Z][A-Z](?=\\.)")) %>%
                mutate(filename = file.path(dir_target_adm1, paste0(GID_0, ".feather"))) %>%
                as_tibble() %>%
                group_by(GID_0) %>%
                # apply a group_walk to save out (keep = F because we want to keep only GID_1)
                group_walk(~ .x %>% dplyr::select(-filename) %>% write_feather(path = .x$filename[1]), .keep = FALSE)
              
              # aggregate to ADM0 level
              df_adm1 %>%
                # add the ensemble variable in which is not included in the results returned by project_impacts_initial_submission()
                mutate(ensemble = ensemble) %>%
                aggregate_adm1_to_adm0(adm1_weights_used = gdp_weights, is_mc_output = F, is_input_datatable = T) %>%
                write_feather(file.path(dir_target_adm0, paste0(identifier, ".feather")))
              
              return(TRUE)
              
            },
            .options = furrr_options(seed = NULL),
            .progress = TRUE)

# aggregate from ADM0 to global
dir_target_adm0 %>% list.files(recursive = T, full.names = T) %>% map_dfr(read_feather) %>%
  aggregate_adm0_to_global(adm0_weights_used = df_ssp, is_mc_output = FALSE) %>%
  write_feather(file.path(dir_output_data, "pointestimates_stagg_oldprojectimpacts", "df_global_bc.feather"))

