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

# load global warming level data for model-scenario-specific baseline periods based on global warming level (defaulting to 0.84 in src/utils/project_impacts())
df_gwl <- read_csv(file.path("data", "input", "df_gwl.csv"))

# GDP weights for ADM1-to-ADM0 aggregation
gdp_weights <- read_csv(file.path("data", "adm1-weights.csv"))  %>%
  group_by(GID_0) %>% mutate(adm1_weight = GDP/sum(GDP)) %>% ungroup() %>%
  select(GID_0, GID_1, GDP, adm1_weight)

# draw coefficients
# NOTE: get_coefficients() uses a default seed for Monte Carlo draws from multi-variate Gaussian, so we do not need to fix the seed here)
coefs_mc <- get_coefficients(n = 1000, draw_from_vcov = T)
coefs_mc_T_Pt_Pt2zero <- get_coefficients(model_used = "T_Pt_Pt2zero", draw_from_vcov = T, n = 1000)

# map the directory where climate projection files are stored
dir_climate_data_adm1 <- file.path("/net", "cfc", "landclim1", "pwaidelich", "cmip6_ng_adm1_stagg")

# map the output directory where GDP projection files (ADM0-level) will be stored
dir_output_data <- file.path("/net", "cfc", "landclim1", "pwaidelich")

# load the tibble objects storing model-scenarios for which we have all 6 climate indicators available
df_modelscen <- read_csv(file.path("data", "input", "df_modelscen.csv"))

# read in the cap for monthly precip deviation from src/00 Create Projections.R
ceiling_Wn_bc <- readRDS(file.path("data", "ceiling_Wn_bc.rds"))

# initiate parallel processing
plan(multisession(workers = 40))



################################################################################
############### CREATE PROJECTIONS #############################################
################################################################################

# map the target directory and create it if necessary
if(!dir.exists( file.path(dir_output_data, "fulldistr_stagg"))) dir.create( file.path(dir_output_data, "fulldistr_stagg"))

dir_target_adm0 <- file.path(dir_output_data, "fulldistr_stagg", "adm0")

if(!dir.exists(dir_target_adm0)) dir.create(dir_target_adm0)

# loop through all model runs in df_modelscen and apply produce_gdp_projections(), which is defined in src/utils/project_impacts.R
future_pmap(.l = df_modelscen,
            .f = ~ produce_gdp_projections(model = ..1,
                                           scenario = ..2,
                                           ensemble = ..3,
                                           dir_climate_adm1 = dir_climate_data_adm1,
                                           dir_out_mc = dir_target_adm0, 
                                           overwrite_previous_files = F,
                                           coefs_matrix_used = coefs_mc,
                                           coefs_matrix_used_Tonly = coefs_mc_T_Pt_Pt2zero,
                                           adm1_weights_used = gdp_weights,
                                           df_gwl_baseline = df_gwl,
                                           use_gwl_baseline = TRUE,
                                           ceiling_precipdeviation = ceiling_Wn_bc,
                                           include_tx5d = FALSE),
            .options = furrr_options(seed = NULL),
            .progress = TRUE)