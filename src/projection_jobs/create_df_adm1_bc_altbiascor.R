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
coefs_point <- get_coefficients()
coefs_point_T_Pt_Pt2zero <- get_coefficients(model_used = "T_Pt_Pt2zero")

# map the directory where climate projection files are stored
dir_climate_data_adm1 <- file.path("/net", "cfc", "landclim1", "pwaidelich", "cmip6_ng_adm1_stagg_altbiascor")

# map the output directory where GDP projection files (ADM0-level) will be stored
dir_output_data <- file.path("/net", "cfc", "landclim1", "pwaidelich")

# load the tibble objects storing model-scenarios for which we have all 6 climate indicators available
df_modelscen <- read_csv(file.path("data", "input", "df_modelscen.csv"))

# read in the cap for monthly precip deviation from src/00 Create Projections.R
ceiling_Wn_bc <- readRDS(file.path("data", "ceiling_Wn_bc.rds"))



################################################################################
############### CREATE PROJECTIONS #############################################
################################################################################

# ensure target directories exist
dir_target <- file.path(dir_output_data, "pointestimates_stagg_altbiascor")

# specify the identifier for the example model run used 
identifier <- "MPI-ESM1-2-LR_ssp370_r1i1p1f1"

# calculate ADM1-level impacts using point estimates
df_adm0 <- collect_climdata_feathers(model = "MPI-ESM1-2-LR",
                                    scenario = "ssp370",
                                    ensemble = "r1i1p1f1",
                                    dir_files = dir_climate_data_adm1) %>%
                # apply the ceiling for bias-corrected monthly precip deviation
                mutate(Wn = if_else(Wn > ceiling_Wn_bc, ceiling_Wn_bc, Wn)) %>%
                # calculate GDP impacts
                project_impacts(coefs = coefs_point[, 1],
                                coefs_Tonly = coefs_point_T_Pt_Pt2zero[, 1],
                                return_helper_columns = T,
                                df_gwl_baseline = df_gwl,
                                use_gwl_baseline = T,
                                include_tx5d = F,
                                return_as_datatable = T
                ) %>%
                # discard helper columns that are not of interest
                dplyr::select(-c(imp_temp_during_base, Tmean_fd, dimp_temp)) %>%
                # aggregate to ADM0 level
                aggregate_adm1_to_adm0(adm1_weights_used = gdp_weights, is_mc_output = F, is_input_datatable = T)

# save out ADM0-level results
df_adm0 %>% write_feather(file.path(dir_target, paste0("adm0_", identifier, ".feather")))

# aggregate from ADM0 to global
df_adm0 %>% aggregate_adm0_to_global(adm0_weights_used = df_ssp, is_mc_output = FALSE) %>%
  write_feather(file.path(dir_target, paste0("global_.feather")))

