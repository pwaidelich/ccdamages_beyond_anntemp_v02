# clean out the environment
rm(list = ls())

# load packages
library(tidyverse)       # for general data wrangling & plotting
library(furrr)           # for parallelizing
library(feather)         # for reading in and writing out large files faster

# load helper functions
source(file.path("src", "utils", "analyse_impacts.R"))

# load global warming level data for model-scenario-specific baseline periods
df_gwl <- read_csv(file.path("data", "input", "df_gwl.csv"))

# map where GDP projection files (ADM0-level) will be stored
dir_output_data <- file.path("/net", "cfc", "landclim1", "pwaidelich")

# load model run overview
df_modelscen <- read_csv(file.path("data", "input", "df_modelscen.csv"))

# load list of sovereign ADM0-level countries considered created in a previous script
sovereign_countries_gid0 <- readRDS(file.path("data", "sovereign_countries_gid0.rds"))

# initiate parallel processing
plan(multisession(workers = 40))



################################################################################
############### POINT ESTIMATE SUMMARY STATS & GCM WEIGHTS #####################
################################################################################

# combine all files using point estimates for dose-response functions into one tibble
df_global_bc_pointdistr_le <- read_feather(file.path(dir_output_data, "pointestimates_stagg", "df_global_bc.feather"))

# check # of model-scenario-ensemble pairings
if(!nrow(df_global_bc_pointdistr_le %>% count(model, scenario, ensemble)) == 199) stop("Expected 199 model-scenario-ensemble pairings but got a different no.")

# save out the entire distribution as a .feather file
write_feather(df_global_bc_pointdistr_le, file.path("data", "df_global_bc_pointdistr_le.feather"))

# merge in warming levels
df_global_bc_pointdistr_le_gwl <- df_global_bc_pointdistr_le %>%
  left_join(df_gwl, by = c("model", "scenario", "ensemble", "year")) %>%
  # discard 0.38 & 0.84 GWL and years outside our GWL windows
  filter(warming_level >= 1) %>%
  # rename CESM2-LE to CESM2 to ensure that we treat the single run of CESM2 for SSP1-2.6 as the same model as the large ensemble
  # NOTE: this matters for calculationg GCM weights
  mutate(model = if_else(model == "CESM2-LE", "CESM2", model))

# ensure that we did not lose any model-scenario-ensemble pairings or GWLs
if(!nrow(df_global_bc_pointdistr_le_gwl %>% count(model, scenario, ensemble)) == 199) stop("Expected 199 model-scenario-ensemble pairings but got a different no.")
if(mean(unique(df_global_bc_pointdistr_le_gwl$warming_level) %in% c(1, 1.5, 2, 3, 4)) != 1) stop("Some GWLs seems to have been list")
   
# save out the entire distribution mapped to GLWs as a .feather file
write_feather(df_global_bc_pointdistr_le_gwl, file.path("data", "df_global_bc_pointdistr_le_gwl.feather"))

# weight each GCM inversely by no. of obs per GWL
df_gcm_weights_le <- df_global_bc_pointdistr_le_gwl %>%
  group_by(model, warming_level) %>%
  summarise(gcm_weight = 1/n()) %>% ungroup()

# inspect the number of observations which should be 20/40/60 (except for large ensembles)
df_gcm_weights_le %>% mutate(n = 1/gcm_weight) %>% count(n, warming_level)
# NOTE: CESM2 SSP1-2.6 reaches +2C, MPI-LR SSP1-2.6 reaches 1.5C, hence we have an additional 20
#       on top of the large ensemble size for all GWLs equal or below to that
# NOTE: there is one realization for MPI LE that does not reach +4C, hence 20 obs less

# save out the weights
saveRDS(df_gcm_weights_le, file.path("data", "df_gcm_weights_le.rds"))

# summarize the distribution using summarystats_for_distribution(), which is defined in src/utils/analyse_impacts.R
df_global_bc_pointdistr_le_aggr <- df_global_bc_pointdistr_le_gwl %>%
  # subset to variables of interest
  select(model, scenario, ensemble, year, warming_level, imp_total, imp_varextremes, ends_with("diff")) %>%
  summarystats_for_distribution(gcm_weights_used = df_gcm_weights_le,
                                is_input_adm0 = F,
                                warming_levels_used = c(1, 1.5, 2, 3, 4))

# save out summary stats
write_feather(df_global_bc_pointdistr_le_aggr, file.path("data", "df_global_bc_pointdistr_le_aggr.feather"))

# # quick visualization for checking - UNCOMMENT IF REQUIRED
# df_global_bc_pointdistr_le_aggr %>% 
#   select(warming_level, ends_with("mean"), ends_with("perc10"), ends_with("perc90")) %>%
#   pivot_longer(cols = starts_with("imp")) %>%
#   mutate(summarystat = str_extract(name, "(?<=_)[:alnum:]+$"),
#          name = str_remove(name, "_[:alnum:]+$")) %>%
#   pivot_wider(names_from = "summarystat", values_from = "value") %>%
#   ggplot(aes(as.character(warming_level))) + geom_point(aes(y=mean)) +
#   geom_errorbar(aes(ymin = perc10, ymax = perc90)) +
#   facet_wrap(~ name, scales = "free_y")

# create a distribution using point estimates that only uses first realization per large ensemble ("_cmip6")
df_global_bc_pointdistr_cmip6_gwl <- df_global_bc_pointdistr_le_gwl %>%
  # discard all MPI-ESM1-2-LR ensemble runs except for r1i1p1f1
  filter(!(model == "MPI-ESM1-2-LR" & scenario == "ssp370" & ensemble != "r1i1p1f1")) %>%
  # discard all CESM2 large ensembles runs (= all obs for SSP3-7.0, so no need for scenario filter)
  filter(! (model == "CESM2" & scenario == "ssp370" & ensemble != "r1i1p1f1"))

# ensure that we have the expected number of model runs
if(!nrow(df_global_bc_pointdistr_cmip6_gwl %>% count(model, scenario, ensemble)) == 71) stop("Expected 71 model-scenario-pairings but got different no.")

# save out the distribution as a .feather
write_feather(df_global_bc_pointdistr_cmip6_gwl, file.path("data", "df_global_bc_pointdistr_cmip6_gwl.feather"))

# calculate the weights for inverse GCM sampling for this distribution
df_gcm_weights_cmip6 <- df_global_bc_pointdistr_cmip6_gwl %>%
  group_by(model, warming_level) %>%
  summarise(gcm_weight = 1/n()) %>% ungroup()

# weights should only be 1/20, 1/40 or 1/60 since there are 3 scenarios maximum per model and no large ensembles
if(mean(df_gcm_weights_cmip6$gcm_weight %in% c(1/20, 1/40, 1/60)) != 1) stop("For CMIP6 single runs & 3 RCP-SSPs, GCM weights should be 1/20, 1/40 or 1/60 but are not. Please inspect")

# save out
saveRDS(df_gcm_weights_cmip6, file.path("data", "df_gcm_weights_cmip6.rds"))

# summarize the distribution that excludes large ensembles
df_global_bc_pointdistr_cmip6_aggr <- df_global_bc_pointdistr_cmip6_gwl %>%
  select(model, scenario, ensemble, year, warming_level, imp_total, imp_varextremes, ends_with("diff")) %>%
  summarystats_for_distribution(gcm_weights_used = df_gcm_weights_cmip6,
                                is_input_adm0 = F,
                                warming_levels_used = c(1, 1.5, 2, 3, 4))

# save out
write_feather(df_global_bc_pointdistr_cmip6_aggr, file.path("data", "df_global_bc_pointdistr_cmip6_aggr.feather"))

# clean up the environment
rm(df_global_bc_pointdistr_cmip6_aggr, df_global_bc_pointdistr_cmip6_gwl,
   df_global_bc_pointdistr_le_aggr, df_global_bc_pointdistr_le_gwl, df_global_bc_pointdistr_le)
gc()



################################################################################
################### SUMMARIZE tx5d FILES #######################################
################################################################################

# combine all files into one tibble and discard the single run for CESM2 SSP3-7.0 (since we use large ensemble)
df_global_bc_pointdistr_le_tx5d_gwl <- read_feather(file.path(dir_output_data, "pointestimates_stagg_tx5d", "df_global_bc.feather")) %>%
  left_join(df_gwl, by = c("model", "scenario", "ensemble", "year")) %>%
  filter(warming_level >= 1) %>%
  # rename CESM2-LE to CESM2 to ensure that we treat the single run for different RCP-SSP as the same model as the LE
  mutate(model = if_else(model == "CESM2-LE", "CESM2", model))

# calculate summary stats using weights excl. large ensembles (since we only have Tx5d results for a single run per model)
df_global_bc_pointdistr_le_tx5d_aggr <- df_global_bc_pointdistr_le_tx5d_gwl %>%
  select(model, scenario, ensemble, year, warming_level, imp_total, imp_varextremes, ends_with("diff")) %>%
  # summarize using the weights that feature only one run per large ensemble
  summarystats_for_distribution(gcm_weights_used = df_gcm_weights_cmip6,
                                is_input_adm0 = F,
                                warming_levels_used = c(1, 1.5, 2, 3, 4))

# export
write_feather(df_global_bc_pointdistr_le_tx5d_aggr, file.path("data", "df_global_bc_pointdistr_le_tx5d_aggr.feather"))

# clean up the environment
rm(df_global_bc_pointdistr_le_tx5d_gwl, df_global_bc_pointdistr_le_tx5d_aggr)
gc()



################################################################################
############### CHECK FOR # OF FILES FOR GLOBAL DISTRIBUTION  ##################
################################################################################

# check no. of global files
if((file.path(dir_output_data, "fulldistr_stagg", "global") %>% list.files(recursive = T) %>% length()) != 199000) stop("adm0 should contain 199,000 files but does not. Please inspect")

# put all files going into the full distribution into one vector
subfolders <- file.path(dir_output_data, "fulldistr_stagg", "global") %>% list.files()

if(length(subfolders) != nrow(df_modelscen)) stop("subfolders should have some # of items as df_modelscen has rows")



################################################################################
############### IMPORT & SUMMARIZE  ############################################
################################################################################

# import all files into one tibble and track time required for computation
start <- Sys.time()
df_global_bc_fulldistr_le_gwl <- future_map_dfr(subfolders, ~ file.path(dir_output_data, "fulldistr_stagg", "global", .x) %>%
                                              list.files(full.names = T) %>% map_dfr(read_feather)) %>%
  # merge in GWLs and subset
  left_join(df_gwl, by = c("model", "scenario", "ensemble", "year")) %>%
  filter(!is.na(warming_level) & warming_level >= 1) %>%
  # rename CESM2-LE to CESM2 to ensure that we treat the single run for different RCP-SSP as the same model
  mutate(model = if_else(model == "CESM2-LE", "CESM2", model))

end <- Sys.time()
print("Time for loading and GWL merging for global distribution:")
print(end - start)

rm(start, end)

# ensure that no model-scenario-ensemble is lost through the last step
if(nrow(df_global_bc_fulldistr_le_gwl %>% count(model, scenario, ensemble)) != 199) stop("Merging in GWLs seems to have discarded at least one model-scenario-ensemble combination. Please inspect")

# save out
write_feather(df_global_bc_fulldistr_le_gwl, file.path("data", "df_global_bc_fulldistr_le_gwl.feather"))

# load the GCM weights
df_gcm_weights_le <- readRDS(file.path("data", "df_gcm_weights_le.rds"))

# calculate summary statistics and track time required for computation
start <- Sys.time()
df_global_bc_fulldistr_le_aggr <- df_global_bc_fulldistr_le_gwl %>%
  data.table::as.data.table() %>%
  summarystats_for_distribution(gcm_weights_used = df_gcm_weights_le,
                                is_input_datatable = T,
                                is_input_adm0 = F, warming_levels_used = c(1, 1.5, 2, 3, 4))

end <- Sys.time()
print("Time required for global summary statistics calculation:")
print(end - start)

# save out
write_feather(df_global_bc_fulldistr_le_aggr, file.path("data", "df_global_bc_fulldistr_le_aggr.feather"))
