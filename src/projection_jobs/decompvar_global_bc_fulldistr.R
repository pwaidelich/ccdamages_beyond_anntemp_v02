# load packages
library(tidyverse)
library(dtplyr)
library(feather)

# load the full distribution at the global level
df_global_bc_fulldistr_le_gwl <- read_feather(file.path("data", "df_global_bc_fulldistr_le_gwl.feather"))

# read in the climate model weights
df_gcm_weights_le <- readRDS(file.path("data", "df_gcm_weights_le.rds"))

# prepare the data for variance decomposition
data <- data.table::as.data.table(df_global_bc_fulldistr_le_gwl) %>% lazy_dt() %>%
                                    # subset to warming levels of interest
                                    filter(warming_level %in% c(1.5, 2, 3)) %>%
                                    # create an identifier for each scenario-model run-year pairing
                                    mutate(scenarioensembleyear = paste0(scenario, ensemble, year)) %>%
                                    select(-c(scenario, ensemble, year)) %>%
                                    # include a GID_0 placeholder, so the code can be easily applied to decompose country-level variances if required
                                    mutate(GID_0 = "Global")

# extract median scenario-ensemble-year for each CMIP6 model
df_global_bc_fulldistr_le_medianyear <- data %>% 
  # convert to LONG format
  pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "imp") %>%
  # calculate average impact for each climate indicator ('variable'), model, GWL, and scenario-realization-year
  # NOTE: here, we take within-model averages, so no need for using model weights
  group_by(GID_0, model, warming_level, scenarioensembleyear, variable) %>% summarise(imp = mean(imp)) %>% ungroup() %>%
  # order scenario-realization-years by mean impact and identify median
  group_by(GID_0, model, warming_level, variable) %>%
  arrange(GID_0, model, warming_level, variable, imp) %>%
  mutate(n = n(),
         ind = 1:n()) %>%
  mutate(is_median_scenarioensembleyear = if_else(n %% 2 == 0, ind == n/2, ind == n/2 + 0.5)) %>%
  # for each model, GWL & variable, extract the median scenario-realization-year
  filter(is_median_scenarioensembleyear) %>%
  ungroup()

# collect results & save out
df_global_bc_fulldistr_le_medianyear %>% as_tibble() %>%
  write_feather(file.path("data", "df_global_bc_fulldistr_le_medianyear.feather"))

# identify the median Monte Carlo draw
df_global_bc_fulldistr_le_medianmc <- data %>%
  # merge in GWL weights (since we average values for the same MC draw, which involves different models)
  left_join(data.table::as.data.table(df_gcm_weights_le), by = c("model", "warming_level")) %>%
  # convert to LONG format
  pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "imp") %>%
  # calculate average impact for GWL, climate indicator and MC draw using climate model weights
  group_by(GID_0, warming_level, variable, monte_carlo_draw) %>% summarise(imp = Hmisc::wtd.mean(imp, weights = gcm_weight, normwt = T)) %>% ungroup() %>%
  # for each GWL and climate indicator, identify the median MC draw
  group_by(GID_0, warming_level, variable) %>%
  arrange(GID_0, warming_level, variable, imp) %>%
  mutate(n = n()) %>%
  mutate(ind = 1:n) %>%
  mutate(is_median_mc = if_else(n %% 2 == 0, ind == n/2, ind == n/2 + 0.5)) %>%
  # extract median MC draws for +3C GWL
  filter(is_median_mc & warming_level == 3) %>%
  ungroup()

# collect results & save out
df_global_bc_fulldistr_le_medianmc %>% as_tibble() %>%
  write_feather(file.path("data", "df_global_bc_fulldistr_le_medianmc.feather"))

# read median scenarioensembleyears and MC draws back in
df_global_bc_fulldistr_le_medianyear <- read_feather(file.path("data", "df_global_bc_fulldistr_le_medianyear.feather"))
df_global_bc_fulldistr_le_medianmc <- read_feather(file.path("data", "df_global_bc_fulldistr_le_medianmc.feather"))

# merge in medians for all models & warming levels (since scenario-years are GWL-specific)
data <- data %>% 
  
  # convert to LONG format
  pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "imp") %>%
  
  # merge in median years
  left_join(data.table::as.data.table(df_global_bc_fulldistr_le_medianyear) %>%
                  dtplyr::lazy_dt() %>% select(GID_0, warming_level, model, scenarioensembleyear, variable, is_median_scenarioensembleyear),
            # NOTE: median years are model- & GWL-specific
            by = c("GID_0", "model", "scenarioensembleyear", "variable", "warming_level")) %>%
  
  # df_global_bc_fulldistr_le_medianyear includes only median years, so merging assigns an NA to all other years - we replace these with FALSE
  mutate(is_median_scenarioensembleyear = replace_na(is_median_scenarioensembleyear, F)) %>%
  
  # merge in the median MC draw and, similarly, replace NAs with FALSE
  left_join(data.table::as.data.table(df_global_bc_fulldistr_le_medianmc) %>%
              dtplyr::lazy_dt() %>%
              select(GID_0, monte_carlo_draw, variable, is_median_mc),
            # NOTE: median MC draw is based on +3C GWL, so neither model- nor GWL-specific
            by = c("GID_0", "monte_carlo_draw", "variable")) %>%
  
  # again, replace NAs (= not median MC draw) with FALSE
  mutate(is_median_mc = replace_na(is_median_mc, F))

# variance between CMIP6 models (fixing MC draw & scenario-ensemble-year at median values)
# NOTE: this is independent of the median model chosen, so we only need to calculate it once
data %>%
  # fix model-realization-year & MC draw at median values by discarding other values
  filter(is_median_scenarioensembleyear & is_median_mc) %>%
  
  # merge in climate model weights to calculate mean/variances across different models
  left_join(data.table::as.data.table(df_gcm_weights_le), by = c("model", "warming_level")) %>%
  
  # calculate marginal variance by climate indicator & GWL
  group_by(GID_0, warming_level, variable) %>%
  summarise(variance = Hmisc::wtd.var(imp, weights = gcm_weight, normwt = T),
            variance_source = "model",
            n = n()) %>% ungroup() %>%
  
  # collect data.table results and export as .feather file
  as_tibble() %>% write_feather(file.path("data", "df_global_bc_fulldistr_le_hsiangvar_model.feather"))

# loop through different choices for the median-like climate model to calculate internal variability & dose-response function uncertainty
# NOTE: we do not parallelize here due to high memory loads
for(median_model_used in c("KACE-1-0-G", "EC-Earth3-Veg", "MPI-ESM1-2-LR", "CESM2")) {
  
  # print out the model name to track progress
  print(median_model_used)
  
  print("Calculating between-Monte-Carlo-draw variance...")
  
  # calculate variance over MC draws (fixing model & scenario-realization-year at median values)
   data %>% filter(model == median_model_used & is_median_scenarioensembleyear) %>%
    
    # calculate marginal variance between MC draws
    # NOTE: we only consider values from a single climate model here, so no need for climate model weights
    group_by(GID_0, warming_level, variable) %>%
    summarise(variance = var(imp),
              variance_source = "monte_carlo_draw",
              n = n()) %>% ungroup() %>%
    
     # collect data.table results and export as .feather file
    as_tibble() %>% write_feather(file.path("data", paste0("df_global_bc_fulldistr_le_hsiangvar_mc_", median_model_used, ".feather")))

  print("Calculating internal variability...")

  # variance between scenario-realization-years (fixing model & MC draw at median values)
  data %>% filter(model == median_model_used & is_median_mc) %>%
    
    # calculate marginal variance between scenario-realization-years
    # NOTE: we only consider values from a single climate model here, so no need for climate model weights
    group_by(GID_0, warming_level, variable) %>%
    summarise(variance = var(imp),
              variance_source = "scenarioensembleyear",
              n = n()) %>% ungroup() %>%
    
    # collect data.table results and export as .feather file
    as_tibble() %>%
    write_feather(file.path("data", paste0("df_global_bc_fulldistr_le_hsiangvar_internal_", median_model_used, ".feather")))

}
