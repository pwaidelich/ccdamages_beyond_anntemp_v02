## RUNTIME: 6min at 40 workers and 500GB RAM

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

# load global warming level data for model-scenario-specific baseline periods
df_gwl <- read_csv(file.path("data", "input", "df_gwl.csv"))

# GDP weights for ADM1-to-ADM0 aggregation
gdp_weights <- read_csv(file.path("data", "adm1-weights.csv"))  %>%
  group_by(GID_0) %>% mutate(adm1_weight = GDP/sum(GDP)) %>% ungroup() %>%
  select(GID_0, GID_1, GDP, adm1_weight)

# # GDP weights for ADM0-to-global aggregation
df_ssp <- readRDS(file.path("data", "df_ssp.rds"))

# draw coefficients
# NOTE: get_coefficients() uses a default seed for Monte Carlo draws from multi-variate Gaussian, so we do not need to fix the seed here)
coefs_point_Tx5d <- get_coefficients(model_used = "Kotz + tx5d")
coefs_point_T_Pt_Pt2zero <- get_coefficients(model_used = "T_Pt_Pt2zero")

# map the directory where climate projection files are stored
dir_climate_data_adm1 <- file.path("/net", "cfc", "landclim1", "pwaidelich", "cmip6_ng_adm1_stagg")

# map output directory where GDP projection files (ADM0-level) will be stored
dir_output_data <- file.path("/net", "cfc", "landclim1", "pwaidelich")

# load the tibble objects storing model-scenarios for which we have all 6 climate indicators available
df_modelscen <- read_csv(file.path("data", "input", "df_modelscen.csv"))

# read in the cap for monthly precip deviation from a previous script
ceiling_Wn_bc <- readRDS(file.path("data", "ceiling_Wn_bc.rds"))

# initiate parallel processing
plan(multisession(workers = 35, gc = TRUE))



################################################################################
############### CREATE PROJECTIONS #############################################
################################################################################

# ensure ADM1- and ADM0-level target directories exist
dir_adm1 <- file.path(dir_output_data, "pointestimates_stagg_tx5d", "adm1")
if(!dir.exists(dir_adm1)) dir.create(dir_adm1) 

dir_target_adm0 <- file.path(dir_output_data, "pointestimates_stagg_tx5d", "adm0")
if(!dir.exists(dir_target_adm0)) dir.create(dir_target_adm0)  

# NOTE: the following model-scenario pairs are not projected out because required tasmax is not available
# 1 ssp126   CESM2
# 2 ssp370   CESM2 (= large ensemble)
# 3 ssp370   CESM2-WACCM
# 4 ssp370   CMCC-CM2-SR5

# loop through the model runs in df_modelscen excl. the ones listed above
future_pmap(.l = df_modelscen %>%
              # discard the models listed above
              filter(!str_detect(model, "CESM2")) %>% 
              filter(!(model == "CMCC-CM2-SR5" & scenario == "ssp370")) %>%
              # also discard MPI large ensemble (Tx5d currently not available)
              filter(!(model == "MPI-ESM1-2-LR" & scenario == "ssp370" & ensemble != "r1i1p1f1"))
              ,
            .f = function(model, scenario, ensemble) {
              
              # create identifier character
              identifier <- paste0(model, "_", scenario, "_", ensemble)
              
              # print out progress - commented out since this only works for map(), not future_map()
              # paste0("Calculating point estimate results for ", identifier) %>% print()
              
              # calculate ADM1-level impacts using point estimates
              df_adm1 <- collect_climdata_feathers(model = model,
                                                   scenario = scenario,
                                                   ensemble = ensemble,
                                                   dir_files = dir_climate_data_adm1,
                                                   include_tx5d = T) %>%
                # apply the ceiling for bias-corrected monthly precip deviation
                mutate(Wn = if_else(Wn > ceiling_Wn_bc, ceiling_Wn_bc, Wn)) %>%
                # calculate GDP impacts
                project_impacts(coefs = coefs_point_Tx5d[, 1],
                                coefs_Tonly = coefs_point_T_Pt_Pt2zero[, 1],
                                return_helper_columns = T,
                                df_gwl_baseline = df_gwl,
                                use_gwl_baseline = T,
                                include_tx5d = T,
                                return_as_datatable = T
                ) %>%
                # discard helper columns that are not of interest
                dplyr::select(-c(imp_temp_during_base, Tmean_fd, dimp_temp))
              
              # determine the target directory for ADM1-level results
              dir_target_adm1 <- file.path(dir_adm1, identifier)
              if(!dir.exists(dir_target_adm1)) dir.create(dir_target_adm1)
              
              # write out in ADM0-level chunks (for easy access)
              df_adm1 %>% as_tibble() %>%
                mutate(GID_0 = str_extract(GID_1, "^[A-Z][A-Z][A-Z](?=\\.)")) %>%
                mutate(filename = file.path(dir_target_adm1, paste0(GID_0, ".feather"))) %>%
                group_by(GID_0) %>%
                group_walk(~ .x %>% dplyr::select(-filename) %>% write_feather(path = .x$filename[1]), .keep = FALSE)
              
              # aggregate to ADM0 level
              df_adm1 %>%
                aggregate_adm1_to_adm0(adm1_weights_used = gdp_weights, is_mc_output = F, is_input_datatable = T) %>%
                write_feather(file.path(dir_target_adm0, paste0(identifier, ".feather")))
              
              return(TRUE)
              
            }
)

# aggregate from ADM0 to global
dir_target_adm0 %>% list.files(recursive = T, full.names = T) %>% map_dfr(read_feather) %>%
  aggregate_adm0_to_global(adm0_weights_used = df_ssp, is_mc_output = FALSE) %>%
  write_feather(file.path(dir_output_data, "pointestimates_stagg_tx5d", "df_global_bc.feather"))



