# clean out the environment
rm(list = ls())

# load packages
library(janitor)
library(PupillometryR)   # for raincloud charts
library(ggpattern)       # for hatching in maps created with geom_sf()
library(tidyverse)       # for general data wrangling & plotting
library(ggpubr)          # for grid charts via ggarrange()
library(sf)              # for GIS analysis
library(readxl)          # for importing Excel files
library(furrr)            # for parallel processing
library(feather)          # for importing .feather files
library(ggcorrplot)      # for correlation matrix plots
library(xtable)          # for exporting tables
library(dtplyr)          # to speed up data wrangling for large tibbles
library(ggridges)

# set the ggplot2 theme
theme_set(theme_classic())

# source utils functions
source(file.path("src", "utils", "analyse_impacts.R"))
source(file.path("src", "utils", "project_impacts.R"))
source(file.path("src", "utils", "extract_climate_indices.R"))

# # map the directory where climate projection files are stored (each vector has N_countries x N_modelscenariopairs length)
dir_imp <- file.path("/net", "cfc", "landclim1", "pwaidelich")

# map the Sharepoint directory where climate projection files are stored
dir_climate_data_adm1 <- file.path("/net", "cfc", "landclim1", "pwaidelich", "cmip6_ng_adm1_stagg")

# read in the GCM weights
df_gcm_weights_cmip6 <- readRDS(file.path("data", "df_gcm_weights_cmip6.rds"))
df_gcm_weights_le <- readRDS(file.path("data", "df_gcm_weights_le.rds"))

# read in the data objects required
df_global_bc_pointdistr_le_gwl <- read_feather(file.path("data", "df_global_bc_pointdistr_le_gwl.feather"))
df_global_bc_pointdistr_le_aggr <- read_feather(file.path("data", "df_global_bc_pointdistr_le_aggr.feather"))
df_global_bc_fulldistr_le_aggr <- read_feather(file.path("data", "df_global_bc_fulldistr_le_aggr.feather"))
df_global_bc_fulldistr_le_gwl <- read_feather(file.path("data", "df_global_bc_fulldistr_le_gwl.feather"))
df_adm0_bc_fulldistr_le_aggr <- read_feather(file.path("data", "df_adm0_bc_fulldistr_le_aggr.feather"))

# read in the world shapefile as well as the list of sovereign countries from src/00 Prepare Analysis.R
world <- readRDS(file.path("data", "world.rds"))
sovereign_countries_gid0 <- readRDS(file.path("data", "sovereign_countries_gid0.rds"))

# World Bank data on income-based country classifications
df_worldbank <- read_excel(file.path("data", "input", "230108 World Bank Country Classification.xlsx"),
                           sheet = "List of economies") %>%
  # discard all aggregates (e.g. EU) that have an NA in the region column
  filter(!is.na(Region))

# load the global warming levels
df_gwl <- read_csv(file.path("data", "input", "df_gwl.csv"))

# load the model run overview
df_modelscen <- read_csv(file.path("data", "input", "df_modelscen.csv"))

# deactivate spherical geometry package use to speed up geom_sf() for creating maps
sf_use_s2(FALSE)

# make a LONG-format version of global aggregate stats
df_global_bc_fulldistr_le_aggr_long <- df_global_bc_fulldistr_le_aggr %>%
  setNames(names(.) %>% str_replace("agreement_meansign", "agreementmeansign")) %>%
  pivot_longer(cols = starts_with("imp"), names_to = "variable") %>%
  # extract the summary stat label in the 'variable' column, then remove it
  mutate(summary_stat = str_extract(variable, "_[A-Za-z0-9]+$") %>% str_remove("^_"),
         variable = str_remove(variable, "_[A-Za-z0-9]+$")) %>%
  # convert to wide format and use chart-compatible labels for the climate indicators
  pivot_wider(names_from = "summary_stat", values_from = "value") %>%
  setNames(names(.) %>% str_replace("agreementmeansign", "agreement_meansign"))


################################################################################
############################## ADM1-level RESULTS ##############################
################################################################################

# NOTE: all charts for ADM1-level results are based only on point estimates since
# they only serve for the SI. If Monte Carlo results were available at ADM1-level,
# this would increase storage requirements for outputs by > 1 order of magnitude
# Figures produced in this part of the script feature in Appendix A of the SI


### a) load in all required ADM1-level data objects

# compile the 71 identifiers for single model runs (= excl. large ensembles)
identifier_wo_ensembles <- df_modelscen %>%
  filter(!(model == "CESM2-LE" & ensemble != "r1i1p1f1")) %>%
  filter(!(model == "MPI-ESM1-2-LR" & ensemble != "r1i1p1f1")) %>%
  mutate(identifier = paste0(model, "_", scenario, "_", ensemble)) %>% 
  pull(identifier)

# ensure that indeed, 71 identifiers are selected
if(length(identifier_wo_ensembles) != 71) stop("Should be 71 model runs")

# read in the ADM1-level results using point estimates for dose-response functions
# NOTE: this object is 9 GB large and reading it in can take several minutes
df_adm1 <- file.path(dir_imp, "pointestimates_stagg", "adm1", identifier_wo_ensembles) %>% list.files(full.names = T) %>%
  map_dfr(read_feather) %>%
  data.table::as.data.table() %>%
  dtplyr::lazy_dt() %>%
  select(-starts_with("imp"), -starts_with("dimp"), -Tmean_base) %>%
  as_tibble()

# load the data from Kotz et al 2022
df_kotz2022 <- read.csv(file.path("..", "sharepoint", "Data", "Kotz et al - data & code", "Data_code_for_Rainfall_economic_production_Kotz_2021", "zenodo", "secondary_data", "PT_master.csv")) %>%
  # create the same ADM1-level identifier
  mutate(GID_1 = paste0(iso, ".", id_1, "_1")) %>%
  # filter to baseline period 1979--2019
  filter(year >= 1979 & year <= 2019) %>% as_tibble()

# filter down to the +0.84 baseline period and merge in Kotz et al. 2022 sample data
df_adm1_comparison <- df_adm1 %>% filter(is_baseline) %>%
  # subset to the regions covered in Kotz et al. 2022
  mutate(in_kotz_sample = GID_1 %in% df_kotz2022$GID_1) %>%
  as_tibble()

# write this object out, so we do not need to reload the full ADM1-level results
write_feather(df_adm1_comparison, file.path("data", "df_adm1_comparison.feather"))


# read back in - UNCOMMENT TO AVOID RE-CREATING df_adm1_comparison ABOVE
# df_adm1_comparison <- read_feather(file.path("data", "df_adm1_comparison.feather"))

# read in all ERA5 files and merge them together
df_era5_adm1 <- map(file.path(dir_climate_data_adm1, "era5_084") %>%
                      list.files(full.names =T) %>% str_subset("mon_totd", negate = T),
                    .f = ~ read_feather(.x)) %>%
  # combine the list produced by map() into a single data frame via a left_join
  reduce(left_join, by = c("year", "GID_1"))

# calculate annual precipitation based on monthly totals
df_era5_Pt <- file.path(dir_climate_data_adm1, "era5_084") %>%
  list.files(full.names =T) %>% str_subset("mon_totd") %>% read_feather() %>%
  mutate(year_clean = lubridate::year(lubridate::ym(year))) %>%
  group_by(GID_1, year_clean) %>% summarise(Pt = sum(mon_totd)) %>%
  rename(year = "year_clean")

# merge the two ERA files together and merge in values from Kotz et al 2022 for benchmarking
df_era5_adm1 <- df_era5_adm1 %>% left_join(df_era5_Pt, by = c("GID_1", "year")) %>%
  left_join(df_kotz2022 %>% dplyr::select(GID_1, year, Tmean, Tstd, Pt, Wn, wet_days_1, vwet_days1_am_99p9),
            by = c("GID_1", "year"), suffix = c("", "_kotz"))


### b) compare baseline period distributions

# plot density for ERA5 & CMIP6 baseline
densities <- map(c("Tmean", "Tstd", "Pt", "Wn", "wet_days_1", "vwet_days1_am_99p9"),
                 .f = function(var_selected) {
                   df_adm1_comparison %>%
                     # NOTE: uncomment the line below to test the chart at lower sample size if required
                     # slice_sample(n = round(nrow(df_adm1_comparison)/50, 0)) %>%
                     pivot_longer(cols = all_of(var_selected),
                                  names_to = "variable", values_to = "value") %>%
                     mutate(variable = rewrite_climate_label_with_units(variable)) %>%
                     ggplot(aes(x = value)) +
                     geom_histogram(aes(y = ..density.., colour = "CMIP6 baseline (+0.84\u00B0C GWL)", fill = "CMIP6 baseline (+0.84\u00B0C)"), alpha = 0.3) +
                     geom_histogram(data = df_era5_adm1  %>% pivot_longer(cols = all_of(var_selected),
                                                                        names_to = "variable", values_to = "value") %>%
                                    mutate(variable = rewrite_climate_label_with_units(variable)),
                                  aes(y = ..density.., x = value, colour ="ERA5 baseline (1979-2019)", fill = "ERA5 baseline (1979-2019)"), alpha = 0.1) + 
                     facet_wrap(~ variable, scales = "free") +
                     theme_classic() +
                     labs(x = NULL, y = "Density", colour = NULL) +
                     guides(fill = "none") +
                     theme(legend.position = "bottom") +
                     theme(axis.title.y = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks.y = element_blank())
                 })

# combine charts into a grid chart and save out - NOTE: this takes some time due to the sample size
density_adm1_vs_era5 <- ggarrange(plotlist = densities, nrow = 2, ncol = 3, common.legend = T, legend = "bottom", align = "hv")  
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " ADM1distribution_cmip6_vs_era5.png")),
       plot = density_adm1_vs_era5,
       width = 6.5, height = 4)

# make ridge density chart for each model separately
figure_ridgedensity_cmip6era5 <- function(var_selected = NULL) {
  df_adm1_comparison %>%
    #slice_sample(n = nrow(df_adm1_comparison)) %>%
    select(model, var_for_plot = var_selected) %>%
    bind_rows(df_era5_adm1 %>% mutate(model = "ERA5") %>% select(model, var_for_plot = var_selected)) %>%
    mutate(model = factor(model, levels = c(unique(df_adm1_comparison$model) %>% sort(), "ERA5"))) %>%
    ggplot(aes(var_for_plot, model, height = after_stat(density))) + ggridges::geom_density_ridges(aes(fill = model=="ERA5"), stat = "density", trim = T) +
    labs(x = rewrite_climate_label_with_units(var_selected), y = NULL) +
    theme(legend.position = "none")
}

# NOTE: this plot takes some time to compute due to the large size of df_adm1_comparison
density_cmip6_vs_era5_bymodel <- ggarrange(plotlist = map(c("Tmean", "Tstd", "Pt", "wet_days_1", "Wn", "vwet_days1_am_99p9"), figure_ridgedensity_cmip6era5),
                                           nrow = 2, ncol = 3, align = "hv")
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " ADM1distribution_cmip6_vs_era5_bymodel.pdf")),
       plot = density_cmip6_vs_era5_bymodel,
       width = 10, height = 10)

# clean up
rm(density_adm1_vs_era5, density_adm1_vs_era5, densities)
gc()

# calculate baseline correlation coefficients for CMIP6
cor_climate_adm1_cmip6 <- df_adm1_comparison %>% select(Tmean, Tstd, Pt, Wn, wet_days_1, vwet_days1_am_99p9) %>%
  setNames(names(.) %>% rewrite_climate_label()) %>%
  cor(use = "pair")

# calculate baseline correlation coefficients for ERA5
cor_climate_adm1_era5 <- df_era5_adm1 %>%  select(Tmean, Tstd, Pt, Wn, wet_days_1, vwet_days1_am_99p9) %>%
  setNames(names(.) %>% rewrite_climate_label()) %>%
  cor(use = "pair")

# plot the two corrplots next to each for comparisons and save out
ggarrange(
  ggcorrplot(cor_climate_adm1_era5, outline.col = "white", type = "lower", hc.order = T,
             ggtheme = theme_classic, lab = T, pch.col = "grey", pch.cex =  10, digits = 2) +
    labs(subtitle = "ERA5 baseline (1979-2019)"),
  ggcorrplot(cor_climate_adm1_cmip6, outline.col = "white", type = "lower", hc.order = T,
             ggtheme = theme_classic, lab = T, pch.col = "grey", pch.cex =  10, digits = 2) +
    labs(subtitle = "CMIP6 baseline (+0.84\u00B0C GWL)"),
  nrow = 1, ncol = 2
)
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " ADM1correlation_era5_vs_cmip6.png")),
       width = 12, height = 10)


### c) benchmark our ERA5-based values against the sample from Kotz et al. (2022)

# baseline distribution of our ERA5 values vs Kotz et al 2022 ERA5 values
df_era5_adm1 %>%
  filter(!is.na(Tmean) & !is.na(Tmean_kotz)) %>%
  select(-contains("tx5d")) %>%
  # NOTE: uncomment the line below to test the chart at lower sample size if required
  #slice_sample(n = round(nrow(df_adm1_comparison)/20, 0)) %>%
  pivot_longer(cols = c(starts_with("T"), starts_with("Pt"), contains("wet_day"), starts_with("Wn")),
                                    names_to = "variable", values_to = "value") %>%
  mutate(is_kotz = if_else(str_detect(variable, "_kotz"),
                           "Kotz et al. (2022), 1979-2019",
                           "ERA5, 1979-2019 (our values)"),
         variable = str_remove(variable, "_kotz$")) %>%
  mutate(variable = rewrite_climate_label_with_units(variable)) %>%
  ggplot(aes(x = value)) + geom_density(aes(colour = is_kotz, fill = is_kotz), alpha = 0.2) +
  facet_wrap(~ variable, scales = "free") +
  theme_classic() +
  labs(x = NULL, y = "Density", colour = NULL) +
  guides(fill = "none") +
  theme(legend.position = "bottom")
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " ADM1distribution_era5_vs_kotz2022.png")),
       width = 6.5, height = 4)

# correlation coefficients between our ERA5 values and Kotz et al. (2022)
tibble(variable = c("Tmean", "Tstd",
                    "Wn", "wet_days_1",
                    "vwet_days1_am_99p9", "Pt"),
       cor = c(
         cor(df_era5_adm1$Tmean, df_era5_adm1$Tmean_kotz, use = "pair"),
         cor(df_era5_adm1$Tstd, df_era5_adm1$Tstd_kotz, use = "pair"),
         cor(df_era5_adm1$Wn, df_era5_adm1$Wn_kotz, use = "pair"),
         cor(df_era5_adm1$wet_days_1, df_era5_adm1$wet_days_1_kotz, use = "pair"),
         cor(df_era5_adm1$vwet_days1_am_99p9, df_era5_adm1$vwet_days1_am_99p9_kotz, use = "pair"),
         cor(df_era5_adm1$Pt, df_era5_adm1$Pt_kotz, use = "pair")
       )) %>%
  mutate(variable = rewrite_climate_label(variable)) %>%
  mutate(cor = round(cor, 3)) %>%
  setNames(names(.) %>% str_replace("^variable$", "Climate indicator") %>%
             str_replace("^cor$", "Correlation coefficients")) %>%
  gt() %>%
  tab_header(
    title = md("**Correlation of ERA5-based climate indicators between our approach and Kotz et al. (2022)**"),
    subtitle = md("(calculated for annual, ADM1-level indicators across 1979-2019)")
  )




################################################################################
######################## OVERVIEW TABLE: GLOBAL IMPACTS ########################
################################################################################

# print out means as well as upper/lower deciles for all impact channels & global warming levels
df_global_bc_fulldistr_le_aggr %>%
  mutate(imp_temp_diff = paste0(100*round(imp_temp_diff_mean, 4), "% (",
                                100*round(imp_temp_diff_perc10, 4), " to ",
                                100*round(imp_temp_diff_perc90, 4), "%)"),
         imp_Pt_diff = paste0(100*round(imp_Pt_diff_mean, 4), "% (",
                              100*round(imp_Pt_diff_perc10, 4), " to ",
                              100*round(imp_Pt_diff_perc90, 4), "%)"),
         imp_Tstd_diff = paste0(100*round(imp_Tstd_diff_mean, 4), "% (",
                                100*round(imp_Tstd_diff_perc10, 4), " to ",
                                100*round(imp_Tstd_diff_perc90, 4), "%)"),
         imp_Wn_diff = paste0(100*round(imp_Wn_diff_mean, 4), "% (",
                              100*round(imp_Wn_diff_perc10, 4), " to ",
                              100*round(imp_Wn_diff_perc90, 4), "%)"),
         imp_wet_days_1_diff = paste0(100*round(imp_wet_days_1_diff_mean, 4), "% (",
                                      100*round(imp_wet_days_1_diff_perc10, 4), " to ",
                                      100*round(imp_wet_days_1_diff_perc90, 4), "%)"),
         imp_vwet_days1_am_99p9_diff = paste0(100*round(imp_vwet_days1_am_99p9_diff_mean, 4), "% (",
                                              100*round(imp_vwet_days1_am_99p9_diff_perc10, 4), " to ",
                                              100*round(imp_vwet_days1_am_99p9_diff_perc90, 4), "%)"),
         imp_total = paste0(100*round(imp_total_mean, 4), "% (",
                            100*round(imp_total_perc10, 4), " to ",
                            100*round(imp_total_perc90, 4), "%)"),
         imp_temp_Tonly_diff = paste0(100*round(imp_temp_Tonly_diff_mean, 4), "% (",
                                              100*round(imp_temp_Tonly_diff_perc10, 4), " to ",
                                              100*round(imp_temp_Tonly_diff_perc90, 4), "%)"),
         # calculate difference between main results and status quo temperature impacts
         diff_statusquo = paste0(100*round(imp_total_mean - imp_temp_Tonly_diff_mean, 4), "%"),
         # write the warming level as "+3C" instead of "3"
         warming_level = paste0("+", warming_level, "\u00B0C")) %>%
  # select variables of interest
  select(warming_level, imp_temp_diff, imp_Pt_diff,imp_Tstd_diff,imp_Wn_diff, imp_wet_days_1_diff,
         imp_vwet_days1_am_99p9_diff, imp_total, imp_temp_Tonly_diff, diff_statusquo) %>%
  # reformat (one row per indicator, one column per GWL)
  pivot_longer(-warming_level) %>% pivot_wider(names_from=warming_level, values_from=value) %>%
  # use clean indicator names
  mutate(name = name %>% rewrite_impact_label()) %>%
  # rename the 'name' column
  setNames(names(.) %>% str_replace("^name$", "Impact channel")) %>%
  # export as .tex table
  xtable(label = "tab:si_globalresults",
        caption = "Global mean GDP impacts by warming level and climate indicator (upper and lower decile in parentheses)") %>%
  print(file = file.path("..", "sharepoint", "Tables", paste0(Sys.Date(), " si_globalresults_overview_bc.tex")),
        size="\\tiny")


################################################################################
############################ FIGURE 2 ##########################################
################################################################################

### Figure 2 Panels a & b
  
# function to plot distribution & summary stats of global GDP impacts by impact variable
figure_impdistr_global <- function(data_aggr = NULL,
                                   data_distr_gwl = NULL,
                                   var_selection = c("imp_Tstd_diff",
                                                    "imp_Wn_diff",
                                                    "imp_wet_days_1_diff",
                                                    "imp_vwet_days1_am_99p9_diff"),
                                   warming_levels_shown = c(1.5, 2, 3)) {
  
  data_aggr %>%
    # replace agreement_meansign with 'agreement_meansign', such that '_' consistently separates variable name and summary stat
    setNames(names(.) %>% str_replace("agreement_meansign", "agreementmeansign")) %>%
    
    # convert to LONG format
    pivot_longer(cols = starts_with("imp"), names_to = "variable") %>%
    
    # extract the summary stat label in the 'variable' column, such that we have separate columns for impact channel & summary stat
    mutate(summary_stat = str_extract(variable, "_[A-Za-z0-9]+$") %>% str_remove("^_"),
           variable = str_remove(variable, "_[A-Za-z0-9]+$")) %>%
    
    # subset to the climate indicators & warming levels selected
    filter(variable %in% var_selection & warming_level %in% warming_levels_shown) %>%
    
    # convert to wide format and use clean labels for the climate indicators
    pivot_wider(names_from = "summary_stat", values_from = "value") %>%
    mutate(variable = rewrite_impact_label_linebreaks(variable)) %>%
    
    ggplot(aes(paste0(warming_level, "\u00B0C"), colour = variable)) +
    
    # distribution as a violin chart - this requires the full data, not the aggregate summary stats
    geom_violin(data = data_distr_gwl %>% filter(warming_level %in% warming_levels_shown) %>%
                            pivot_longer(cols = starts_with("imp"), names_to = "variable") %>%
                            filter(variable %in% var_selection) %>%
                            mutate(variable = rewrite_impact_label_linebreaks(variable)),
                aes(y = value, colour = variable, fill = variable),
                position = position_dodge(width = 0.6),
                alpha = 0.3,
                linewidth = 0.1) +
    
    # add horizontal dashed line to mark zero
    geom_hline(yintercept = 0, colour = "black", linetype = "dashed", alpha = 0.8)  +
    
    # point for the average impact, errorbar indicating upper-to-lower-decile range
    geom_point(aes(y = mean), position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = perc10, ymax = perc90), position = position_dodge(width = 0.6), width =0.3) +
    
    # suppress certain legend items & axis titles, label axes and set legend style
    guides(fill = NULL, alpha = NULL) + 
    labs(y = "Global GDP impact", x = "Global warming level",colour = NULL, fill = NULL ) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.background = element_blank(),
          legend.box.background = element_blank())
}

# panel for total impacts, annual temp, annual precip and variability & extremes (joint impacts) - with manual GDP shock comparisons
# NOTE: selected year-to-year GDP contractions are sourced from 'data/input/230503 World Bank GDP_2015USD.csv'
# NOTE: plotting the full distribution takes a few minutes
figure2a <- figure_impdistr_global(df_global_bc_fulldistr_le_aggr, df_global_bc_fulldistr_le_gwl,
                                    var_selection = c("imp_total", "imp_temp_diff", "imp_varextremes", "imp_Pt_diff")) +
  geom_hline(yintercept = -0.263, linetype = "dotted", colour = "grey") +
  annotate("text", x = 0.5, y = -0.253, label = "Syrian GDP contraction in 2012 (civil war)", size = 3, fontface = 2, hjust = 0, colour = "azure4") +
  geom_hline(yintercept = -0.1015, linetype = "dotted", colour = "grey") +
  annotate("text", x = 0.5, y = -0.0915, label = "Greece GDP contraction in 2011 (Euro crisis)", size = 3, fontface = 2, hjust = 0, colour = "azure4") +
  scale_colour_manual(values = c("#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#229954")) +
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#229954"))

# same chart for the four variability & extremes indicators disaggregated - again, with manual points of comparison based on World Bank WDI data (for 'World')
figure2b <- figure_impdistr_global(df_global_bc_fulldistr_le_aggr, df_global_bc_fulldistr_le_gwl) +
  geom_hline(yintercept = -0.0134, linetype = "dotted", colour = "grey") +
  annotate("text", x = 0.5, y = -0.0114, label = "Global GDP contraction in 2009 (financial crisis)", size = 3, fontface = 2, hjust = 0, colour = "azure4") +
  geom_hline(yintercept = -0.0311, linetype = "dotted", colour = "grey") +
  annotate("text", x = 0.5, y = -0.0291, label = "Global GDP contraction in 2020 (Covid pandemic)", size = 3, fontface = 2, hjust = 0, colour = "azure4") +
  scale_colour_manual(values = c("#229954", "#229954", "#229954", "#229954"))


### Figure 2 Panels c & d

# write a function to decompose variance following the Schwarzwald & Lenssen (2022) approach
decompose_variance_schwarzwald <- function(data_pointdistr_gwl = NULL, # distribution based on dose-response function point estimates
                                           data_fulldistr_gwl = NULL, # full distribution incl. MC draws for dose-response functions
                                           gcm_weights_used = NULL, # climate model weights
                                           var_rescale_factor = 10^4, # arbitrary factor to rescale impacts to prevent variances from being values like 0.00001%
                                           warming_levels_used = c(1.5, 2, 3) # GLWs displayed
                                           ) {
  
  # if we have no GID_0 information (= results are global), add this variable
  # NOTE: this ensures that the decomposition function works for global & national results
  if(!"GID_0" %in% names(data_pointdistr_gwl) & !"GID_0" %in% names(data_fulldistr_gwl)) {
    data_pointdistr_gwl <- data_pointdistr_gwl %>% mutate(GID_0 = "Global")
    data_fulldistr_gwl <- data_fulldistr_gwl %>% mutate(GID_0 = "Global")
  }
  
  # calculate mean and variance of respective impact by CMIP6 model
  df_within_model <- data.table::as.data.table(data_pointdistr_gwl) %>%
    dtplyr::lazy_dt() %>%
    # select variables required
    select(model, scenario, ensemble, GID_0, year, warming_level, imp_total, imp_varextremes, ends_with("diff"),
           -contains("Tonly")) %>%
    # filter to GWLs of interest
    filter(warming_level %in% warming_levels_used) %>%
    # convert to LONG format and rescale impacts
    pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "imp") %>%
    mutate(imp = imp * var_rescale_factor) %>%
    # calculate mean and variance by climate model for each GWL, territory & impact channel
    group_by(model, warming_level, GID_0, variable) %>%
    summarise(var_within_model = var(imp),
              mean_within_model = mean(imp)) %>% ungroup() %>%
    # merge in climate model weights
    left_join(data.table::as.data.table(gcm_weights_used), by = c("model", "warming_level"))
  
  # calculate variance of means by MC draw using GCM weights
  # NOTE: steps are mostly identical to the ones above
  df_within_mc_draw <- data.table::as.data.table(data_fulldistr_gwl) %>%
    dtplyr::lazy_dt() %>%
    filter(warming_level %in% warming_levels_used) %>%
    pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "imp") %>%
    mutate(imp = imp * var_rescale_factor) %>%
    # merge in climate model weights because we average across different models
    left_join(data.table::as.data.table(gcm_weights_used), by = c("model", "warming_level")) %>%
    # calculate mean impact by MC draw for each GWL, territory and impact channel
    group_by(warming_level, variable, monte_carlo_draw, GID_0) %>%
    summarise(mean_within_mcdraw = Hmisc::wtd.mean(imp, weights = gcm_weight)) %>% ungroup() %>%
    # calculate variance of this mean across MC draws (= dose-response function uncertainty)
    group_by(warming_level, variable, GID_0) %>%
    summarise(damage_function_uncertainty = var(mean_within_mcdraw)) %>% ungroup()
  
  # calculate internal variability & model uncertainty, then merge in dose-response function uncertainty
  df_out <- df_within_model %>% as_tibble() %>%
    # calculate mean of within-model variance (= internal variability) and variance of mean impacts by model (= climate model uncertainty)
    group_by(warming_level, GID_0, variable) %>% 
    summarise(internal_variability = Hmisc::wtd.mean(var_within_model, weights = gcm_weight, normwt = T),
              model_uncertainty = Hmisc::wtd.var(mean_within_model, weights = gcm_weight, normwt = T)) %>% ungroup() %>%
    # merge in dose-response function uncertainty calculated in the previous step
    left_join(df_within_mc_draw %>% as_tibble(), by = c("warming_level", "variable", "GID_0"))
  
  return(df_out)
}

# function to make a figure based on the results returned by decompose_variance_schwarzwald()
figure_vardisaggr_global_schwarzwald <- function(data = NULL, # data object returned by decompose_variance_schwarzwald()
                                                 data_fulldistr_aggr = NULL, # aggregate summary stats for the entire distribution
                                                 var_selection = c("imp_total", "imp_temp_diff", "imp_Pt_diff", "imp_varextremes"), # variables displayed
                                                 use_variance_share = T, # whether variances shares (TRUE) or absolute variance components (FALSE) are displayed
                                                 internal_variability_label = "Internal\nvariability", # label used for within-model variance
                                                 warming_levels_shown = c(1.5, 2, 3), # GWLs displayed
                                                 nrow_used = length(var_selection)/4, # rows used in final chart
                                                 ncol_used = NULL, # cols used in final charts
                                                 var_rescale_factor = 10^4 # variance rescale factor (see above)
                                                 ) {
  
  # if we do not use variance shares, we need to extract total variance from data_fulldistr_aggr and also rescale it
  if(!use_variance_share) {
    
    data_totalvar <- data_fulldistr_aggr %>%
        select(warming_level, ends_with("var")) %>%
        pivot_longer(cols = ends_with("var"), names_to = "variable", values_to = "var") %>%
        mutate(variable = str_remove(variable, "_var$")) %>%
        filter(variable %in% var_selection) %>%
        # we rescale variance with the rescale factor squared since var_rescale_factor is multiplied by impacts above
        mutate(var = var * (var_rescale_factor)^2,
               variable = variable %>% str_remove("_var$") %>% rewrite_impact_label_linebreaks())
  }
  
  # make the chart
  data %>%
    
    # filter to variables of interest
    filter(variable %in% var_selection) %>%
    
    # convert to LONG format
    pivot_longer(cols = c("internal_variability", "model_uncertainty", "damage_function_uncertainty"),
                 names_to = "variance_source", values_to = "var") %>%
    
    # filter to GWLS of interest
    filter(warming_level %in% warming_levels_shown) %>%
    
    # clean up labels
    mutate(variable = rewrite_impact_label_linebreaks(variable),
           variance_source = factor(case_when(variance_source == "internal_variability" ~ internal_variability_label,
                                              variance_source == "model_uncertainty" ~ "Climate model\nuncertainty",
                                              variance_source == "damage_function_uncertainty" ~ "Dose-response function\nuncertainty",
                                              TRUE ~ NA_character_),
                                    levels = c("Dose-response function\nuncertainty", internal_variability_label, "Climate model\nuncertainty"))) %>%
    
    # convert to shares if use_variance_share is TRUE
    {if(use_variance_share) group_by(., variable, warming_level) %>% mutate(var = var/sum(var)) %>% ungroup() else .} %>%
    
    # make the bar chart
    ggplot(aes(paste0(warming_level, "\u00B0C"), var)) +
    geom_col(aes(fill = variance_source), alpha = 0.7) +
    
    # include total variance as a triangle point if use_variance_share is FALSE
    {if(!use_variance_share) geom_point(data = data_totalvar %>% filter(warming_level %in% warming_levels_shown), aes(shape = "Total variance\nincl. interaction"), colour = "darkgrey", size = 3)} +
    
    # otherwise label the shares
    {if(use_variance_share) geom_label(aes(fill = variance_source, label = paste0(100*(round(var, ifelse(round(var, 2) == 0, 4, 2))), "%")),
                                       position = position_stack(vjust = 0.5, reverse = F),
                                       size = 2.5,
                                       show.legend = F)} +
    
    # one facet per climate indicator
    facet_wrap(~ variable, scales = ifelse(use_variance_share, "fixed", "free_y"), nrow = nrow_used, ncol = ncol_used) +
    
    # set labels, legends and scales
    labs(y = ifelse(use_variance_share, "Share in variance", "Variance of GDP impacts (in squared bps)"), x = "Global warming level",
         colour = NULL, fill = NULL, shape = NULL) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order=2), shape = guide_legend(order=3)) +
    {if(use_variance_share) scale_y_continuous(labels = scales::percent)} +
    scale_fill_manual(values = c("#F8766D", "#619CFF", "#00BA38")) +
    scale_shape_manual(values = c(17)) +
    theme_classic() +
    theme(legend.position = "bottom")
}

# decompose the variance for the full distribution incl. LEs using Schwarzwald & Lenssen approach
df_var_bc_fulldistr_le_schwarzwald <- decompose_variance_schwarzwald(df_global_bc_pointdistr_le_gwl, df_global_bc_fulldistr_le_gwl,
                                                                     gcm_weights_used = df_gcm_weights_le)

# plot relative variance shares excl. interaction term
figure2c <- df_var_bc_fulldistr_le_schwarzwald %>% figure_vardisaggr_global_schwarzwald()
figure2d <- df_var_bc_fulldistr_le_schwarzwald %>%
  figure_vardisaggr_global_schwarzwald(var_selection = c("imp_Tstd_diff", "imp_Wn_diff", "imp_wet_days_1_diff", "imp_vwet_days1_am_99p9_diff"))


### Figure 2 - panels combined

figure2 <- ggarrange(ggarrange(figure2a, figure2c, nrow = 2, ncol = 1,
                               labels = c("a", "c")),
                     ggarrange(figure2b, figure2d, nrow = 2, ncol = 1,
                               labels = c("b", "d")),
                     nrow = 1, ncol = 2)
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " Figure2.pdf")),
       plot = figure2,
       width = 11, height = 7)

# clean out the environment
rm(figure2a, figure2b, figure2c, figure2d, figure2)
gc()



################################################################################
########################### FIGURE 3 ###########################################
################################################################################

# plot total impacts as a map
figure_map_impacts <- function(data = NULL, # summary stats data
                               var_selected = "imp_total", # variable to be shown in map
                               warming_level_selected = 3, # GWL used for map
                               subtitle_selected = NULL) {
  
  # take world shapefile and remove Antarctica for visual purposes (not included in our sample)
  world %>% filter(sovereignt != "Antarctica") %>% 
    
    # merge in data of interest
    left_join(data %>% filter(warming_level == warming_level_selected,
                              variable == var_selected
                              ), by = "GID_0") %>%
    
    # convert to shapefile
    st_as_sf() %>%
    
    # make the map
    ggplot() + geom_sf(aes(fill = ifelse(GID_0 %in% sovereign_countries_gid0, mean, NA)), size = 0.1) +
    
    # set the colour scale
    scale_fill_gradient2(low = "red", mid = "gray88", high = "blue",
                         n.breaks = 4,
                         labels = scales::percent,
                         # set the upper bound of the fill scale legend to zero if all values in the data are negative
                         limits = c(NA, ifelse((data %>% filter(warming_level == warming_level_selected,
                                                                variable == var_selected,
                                                                GID_0 %in% sovereign_countries_gid0) %>%
                                                  pull(mean) %>% max()) < 0, 0, NA))) +
    
    # set legend items and axis style
    labs(fill = NULL, subtitle = subtitle_selected) +
    theme_classic() +
    theme(legend.position ="bottom",
          legend.direction = "horizontal",
          axis.line = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank())
} 

# make map for total GDP impacts and impacts of variability and extremes at +3C GWL
figure3a <- figure_map_impacts(df_adm0_bc_fulldistr_le_aggr, "imp_total",
                               subtitle_selected = "Total GDP impacts (all indicators) at +3\u00B0C")
figure3b <- figure_map_impacts(df_adm0_bc_fulldistr_le_aggr, "imp_varextremes",
                               subtitle_selected = "GDP impacts of variability & extremes at +3\u00B0C")

# scatter plot using share of distribution agreeing with the mean's sign
figure_agreementshare_adm0 <- function(data = NULL, data_global_long = NULL, warming_level_selected = 3) {
  
  data %>%
    # merge in income group info from df_worldbank as well as the global results
    left_join(df_worldbank %>% rename(GID_0 = "Code"), by = "GID_0") %>%
    
    # NOTE: we add 'Global' as placeholder value for columns not originally included in data_global_long
    bind_rows(data_global_long %>% mutate(GID_0 = "Global", Region = "Global",  `Income group` = "Global")) %>%
    
    # select the six individual climate indicators (= excluding total impacts and total variability & extremes)
    filter(variable %in% c("imp_temp_diff", "imp_Pt_diff", "imp_Tstd_diff", "imp_Wn_diff",
                           "imp_wet_days_1_diff", "imp_vwet_days1_am_99p9_diff")) %>% 
    
    # clean up labels and make a dummy identifying the global values
    mutate(variable = rewrite_impact_label(variable),
           is_global = GID_0 == "Global") %>%
    
    # subset to GWLs of interest and sovereign countries
    filter(warming_level == warming_level_selected, GID_0 %in% c("Global", sovereign_countries_gid0)) %>%
    
    # Income group in the World Bank file is missing for Venezuela, which is either upper or lower middle income (source: https://publications.iadb.org/en/venezuela-still-upper-middle-income-country-estimating-gni-capita-2015-2021#:~:text=In%20the%202022%20World%20Bank,income%20(GNI)%20of%202013.)
    # therefore, we impute this value, also for Kosovo (based on group of all neighboring countries being Upper middle income)
    mutate(`Income group` = case_when(GID_0 %in% c("VEN", "XKO")  & is.na(`Income group`) ~ "Middle income",
                                      TRUE ~ `Income group`)) %>%
    
    # collapse lower and upper middle income into one group
    mutate(`Income group` = if_else(str_detect(`Income group`, "middle income"), "Middle income", `Income group`)) %>%
    
    # make the scatter plot
    ggplot(aes(mean, agreement_meansign)) +
    
    # add vertical dashed line at 0%GDP impacts
    geom_vline(xintercept = c(0), linetype = "dashed") +
    
    # add IPCC likelihood thresholds
    geom_hline(yintercept = c(0.9, 0.666), linetype = "dotted", colour = "grey") +
    
    # add points
    geom_point(aes(colour = `Income group`, shape = `Income group`, size = `Income group`), alpha = 0.6) +
    
    # one facet per impact channel
    facet_wrap(~ variable, scales = "free_x", nrow = 3, ncol =  2) +
    
    # set axes and styles manually
    scale_x_continuous(labels = scales::percent,
                       breaks = scales::pretty_breaks(3)) +
    scale_colour_manual(values = c("black", scales::hue_pal()(3))) +
    scale_shape_manual(values = c(18, 16, 16, 16)) +
    scale_size_manual(values = c(3, 1.5, 1.5, 1.5)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Mean GDP impact", y = "Confidence in sign of mean GDP impact",
         colour = NULL, shape = NULL, size = NULL) +
    guides(shape = guide_legend(nrow = 2), size = guide_legend(nrow = 2), colour =  guide_legend(nrow = 2)) +
    theme_classic() +
    theme(legend.position = "bottom")
}

# make the chart for +3C GWL
figure3c <- figure_agreementshare_adm0(df_adm0_bc_fulldistr_le_aggr, df_global_bc_fulldistr_le_aggr_long, warming_level_selected = 3)

### combine the 3 panels of Figure 3

figure3 <- ggarrange(ggarrange(ggarrange(figure3a + theme(legend.position = "none"), NULL,
                                         as_ggplot(get_legend(figure3a)), nrow = 3, heights = c(6, -1.2, 1)),
                               ggarrange(figure3b + theme(legend.position = "none"), NULL,
                                         as_ggplot(get_legend(figure3b)), nrow = 3, heights = c(6, -1.2, 1)), nrow = 2, ncol = 1, align = "hv", labels = "auto"),
                     ggarrange(figure3c, labels = "c"),
                     nrow = 1, ncol = 2, widths = c(1.3, 1))
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " Figure3.pdf")),
       plot = figure3,
       width = 8.75, height = 5.75)

# clean out the environment
rm(figure3a, figure3b, figure3c, figure3)


################################################################################
################################## FIGURE 4 ####################################
################################################################################

### Figure 4a: global GDP impacts for our main results & status quo

# compare mean temperature damages for main spec & status quo
figure_global_kotzvsstatusquo <- function(data = NULL, warming_levels_shown = c(1.5, 2, 3)) {
  
  # bind together the results for the main results (= imp_total) & the status quo specification (= imp_temp_Tonly)
  bind_rows(data %>% select(warming_level, imp_total_mean, imp_total_perc10, imp_total_perc90) %>%
              mutate(spec = "All climate indicators (= Kotz et al., 2022)"),
            # NOTE: we rename the status quo damages to imp_total and instead identify them via the new 'spec' column
            data %>% select(warming_level, imp_total_mean = "imp_temp_Tonly_diff_mean",
                            imp_total_perc10 = "imp_temp_Tonly_diff_perc10", imp_total_perc90 = "imp_temp_Tonly_diff_perc90") %>%
              mutate(spec = "Annual temperature only (= status quo)")) %>%
    
    # subset to GWLs of interest and convert spec to a factor to control item order
    filter(warming_level %in% warming_levels_shown) %>%
    mutate(spec = factor(spec, levels = c("Annual temperature only (= status quo)", "All climate indicators (= Kotz et al., 2022)"))) %>%
    
    # create labelled bar chart with total impacts at each warming level & errobar for upper-to-lower-decile range
    ggplot(aes(paste0(warming_level, "\u00B0C"), imp_total_mean, group = spec)) +
    
    # bar for means
    geom_col(aes(fill = spec), position = position_dodge()) +
    
    # text labels for means (rounded to 1st digit of percent)
    geom_text(aes(y = 0.005, label = paste0(100*round(imp_total_mean, 3), "%")), position = position_dodge(width = 1), size = 3) +
    
    # errorbar for upper-to-lower-decile range
    geom_errorbar(aes(ymin = imp_total_perc10, ymax = imp_total_perc90), position = position_dodge(width = 1), width = 0.5) +
    
    # set scales, labels and legend style
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
    labs(y = "Global GDP impact", x = "Global warming level", fill = "Impact projections based on...") +
    theme_classic() +
    theme(legend.position = c(0.385, 0.14), legend.background = element_rect(fill='transparent'))
}

# compare total impacts and status quo impacts for 1.5C, 2C and 3C
figure4a <- figure_global_kotzvsstatusquo(df_global_bc_fulldistr_le_aggr)


## Figure 4b: marginal effects of annual mean temperature for different controls

# write a function to calculate marginal effects of annual mean temp. following Kalkuhl & Wenz 2020 (formula in their Footnote 6)
marginal_effect_temperature <- function(T_0 = NULL, T_shock = 1, coef_used = NULL) {
  out <- coef_used["TmeanD"] + coef_used["TmeanLD"] + T_0*(coef_used["TmeanD.Tmean"] + coef_used["TmeanLD.TmeanL"])
  names(out) <- NULL
  return(out)
}

# vectorize the function across initial temperature ('T_0') 
marginal_effect_temperature <- Vectorize(marginal_effect_temperature, vectorize.args = "T_0")

# write a wrapper for Monte Carlo analysis
marginal_effect_temperature_mc <- function(T_0_range = seq(-10, 30, by = 0.05), T_shock = 1, coef_mc_matrix = NULL) {
  map_dfr(1:ncol(coef_mc_matrix), ~ tibble(T_0 = T_0_range,
                                           damages = marginal_effect_temperature(T_0 = T_0_range, T_shock = T_shock, coef_used = coef_mc_matrix[, .x]),
                                           monte_carlo_draw = .x)
  )
}

# set whether or not to (re-)calculate the uncertainty ranges of marginal effects
# NOTE: since this object is computation-intensive, we save it out to avoid re-runs
redo_marginal_effect_object <- FALSE

# use the function to calculate uncertainty ranges by drawing coefficients from multivariate Gaussian for 1,000 MC draws
if(redo_marginal_effect_object) {
  
  # use the MC wrapper for the main spec by Kotz et al (2022) as well as for the status quo spec ('T_Pt_Pt2zero')
  df_margeffect_uncertainty <- bind_rows(marginal_effect_temperature_mc(coef_mc_matrix = get_coefficients(draw_from_vcov = T,
                                                                                                          model_used = "T_Pt_Tstd",
                                                                                                          n=10^5)) %>% mutate(spec = "T_Pt_Tstd"),
                                         marginal_effect_temperature_mc(coef_mc_matrix = get_coefficients(draw_from_vcov = T,
                                                                                                          model_used = "main",
                                                                                                          n=10^5)) %>% mutate(spec = "main"),
                                         marginal_effect_temperature_mc(coef_mc_matrix = get_coefficients(draw_from_vcov = T,
                                                                                                          model_used = "T_Pt_Pt2zero",
                                                                                                          n=10^5)) %>% mutate(spec = "T_Pt_Pt2zero")
  )
  
  # save out
  write_feather(df_margeffect_uncertainty, file.path("data", "df_margeffect_uncertainty.feather"))

} else {
  
  df_margeffect_uncertainty <- read_feather(file.path("data", "df_margeffect_uncertainty.feather"))

}

# calculate mean marginal effect and 95% CI for each initial temperature step
df_figure4b <- read_feather(file.path("data", "df_margeffect_uncertainty.feather")) %>%
  # use data.table to speed up calculations
  data.table::as.data.table() %>% dtplyr::lazy_dt() %>%
  # calculate mean & 95% CI for each initial temperature step & model specification
  group_by(T_0, spec) %>%
  summarise(mean_margeffect = mean(damages) %>% log_to_gdp_impacts(),
            perc025_margeffect = quantile(damages, 0.025) %>% log_to_gdp_impacts(),
            perc975_margeffect = quantile(damages, 0.975) %>% log_to_gdp_impacts()) %>%
  # collect results as a tibble
  as_tibble() %>%
  # clean up the labels
  mutate(spec = factor(spec %>% str_replace("^T_Pt_Pt2zero$", "Annual precipitation only\n (= status quo)") %>%
                         str_replace("^T_Pt_Tstd$", "+ Temperature variability") %>%
                         str_replace("^main$", "+ Other precipitation indicators\n (= Kotz et al., 2022)"),
                       levels = c("Annual precipitation only\n (= status quo)",
                                  "+ Temperature variability",
                                  "+ Other precipitation indicators\n (= Kotz et al., 2022)"
                       )
  ))


# make a line chart of marginal effects with 95% CI ribbon around them
figure4b <- df_figure4b %>% 
  ggplot(aes(x = T_0, colour = spec)) +
  
  # dashed line at 0%GDP impact
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
  
  # line chart for the mean
  geom_line(aes(y = mean_margeffect, linetype = spec)) +
  
  # 95% CI as a geom_ribbon
  geom_ribbon(data = df_figure4b %>% filter(!str_detect(spec, "Temperature variability")),
              aes(ymin = perc025_margeffect, ymax = perc975_margeffect,
                  fill = spec), alpha = 0.15, colour = NA) +
  
  # set scales, legends, and labels
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_colour_manual(values = c("#00BFC4", "black", "#F8766D")) +
  scale_linetype_manual(values = c("solid", "dotted", "solid")) +
  guides(fill = "none") +
  labs(colour = "Controlling for ...", y = "Marginal GDP impact of annual temperature\n(+1\u00B0C increase)",
       fill = "Controlling for ...", linetype = "Controlling for ...",
       x = "Annual temperature (\u00B0C)") +
  theme_classic() +
  theme(legend.position = c(0.3, 0.2), legend.background=element_blank())


### Figure 4c: change of total impacts at +3C vis-a-vis status quo by sovereign ADM0 country

# make a map to compare impacts for all climate indicators and status quo
map_change_vs_statusquo <- function(data = NULL, # ADM0-level summary stats across the full distribution
                                    warming_level_selected = 3, # GWL displayed on the map
                                    subtitle_selected = NULL) {
  
  # take 'world' shapefile, remove Antarctica for visual reasons (not included in our sample)
  df_impactdiff_map <- world %>% filter(!sovereignt == "Antarctica") %>%
    
    # merge in ADM0-level summary stats
    left_join(data %>% filter(warming_level == warming_level_selected,
                              variable %in% c("imp_total", "imp_temp_Tonly_diff")) %>%
                              select(warming_level, GID_0, variable, mean) %>%
                                pivot_wider(names_from = "variable", values_from = "mean") %>%
              # calculate differences between total impacts and status quo impacts
              mutate(diff_statusquo = imp_total - imp_temp_Tonly_diff),
              by = "GID_0") %>%
    
    # set diff_statusquo to NA if ADM0 territory is not sovereign country
    mutate(diff_statusquo = ifelse(GID_0 %in% sovereign_countries_gid0, diff_statusquo, NA))
  
  # make the map
  df_impactdiff_map %>%
    ggplot() +
    geom_sf(aes(fill = diff_statusquo), size = 0.1) +
    scale_fill_gradient2(midpoint = 0, 
                         breaks = c(-0.02, 0),
                         labels = scales::percent,
                         low = "red",
                         # expand the fill scale to zero if necessary
                         limits = c(NA, ifelse(max(df_impactdiff_map$diff_statusquo, na.rm = T) < 0, 0, NA)),
                         high = "blue",
                         mid = "gray88") +
    labs(fill = paste0("Difference in %-pts"), subtitle = subtitle_selected) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          axis.line = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank())
}

# make the maps for difference vs status quo impacts at +3C GWL
figure4c <- map_change_vs_statusquo(df_adm0_bc_fulldistr_le_aggr,
                                    subtitle_selected = "   Change in GDP impact by using all climate indicators (instead of status quo) at +3\u00B0C global warming")


### Figure 4: grid figure on Kotz et al 2022 vs status quo

figure4 <- ggarrange(ggarrange(figure4a, figure4b + labs(y = "Marginal effect of annual temperature \n(+1\u00B0C increase) on GDP"), nrow = 1, ncol = 2, align = "hv", labels = c("a", "b")),
                     ggarrange(figure4c, labels = "c"),
                     heights = c(1, 1.2),
                     nrow = 2, ncol = 1)
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " Figure4.pdf")),
       plot = figure4,
       width = 9, height = 7.5)

# clean out the environment
rm(figure4a, figure4b, figure4c, figure4, df_margeffect_uncertainty)
gc()



################################################################################
########################### FIGURE 5 ###########################################
################################################################################

# load World Bank WDI population figures
df_pop <- read_excel(file.path("data", "input", "231030 World Bank WDI Total Population.xlsx"),
                     sheet = "Data_incl_TWN",
                     na = c("", "NA", ".."))

# calculate the total population for all ADM0 countries included in the analysis
total_pop <- df_pop %>% filter(GID_0 %in% sovereign_countries_gid0 & GID_0 %in% df_adm0_bc_fulldistr_le_aggr$GID_0) %>%
  pull(pop_2022) %>% sum()
total_pop

# calculate the share of world population covered in the analysis
total_pop/(df_pop %>% filter(GID_0 == "WLD") %>% pull(pop_2022))
# -> almost 100%

# merge population info into the tibble with summary statistics at ADM0 country level
df_adm0_with_pop <- df_adm0_bc_fulldistr_le_aggr %>%
  
  # subset to sovereign countries
  filter(GID_0 %in% sovereign_countries_gid0) %>%
  
  # select variables of interest
  select(warming_level, GID_0, variable, ends_with("perc05")) %>%
  
  # merge in population data and calculate weights
  # NOTE: total_pop is not a column, so no need for grouping here
  left_join(df_pop, by = "GID_0") %>%
  mutate(pop_weight = pop_2022/total_pop)

# calculate population share falling below threshold via population-weighted logical
df_pop_exposure <- map_dfr(round(seq(-0.2, 0, by = 0.01), 2),
        .f = function(threshold) {
          
          # for each GWL and impact channel, calculate pop-weighted mean of countries with 5th percentile of impacts under the threshold
          df_adm0_with_pop %>%
            group_by(warming_level, variable) %>% summarise(pop_under_threshold = Hmisc::wtd.mean(perc05 < threshold, weights = pop_weight)) %>%
            # add the threshold used and a spec column
            mutate(perc05_threshold = threshold,
                   spec = if_else(variable == "imp_temp_Tonly_diff", "status_quo", "kotz")) %>%
            ungroup()
        }
)

# write a function to make a line chart of global population shares of exposure to 5\% risks (= using 5th percentile of impacts)
figure_exposure_linechart <- function(data_exposure = NULL, # data on population exposure calculated above
                                      warming_levels_selected = c(1.5, 2, 3) # GWLs of interest
                                      ) {
  
  # subset the exposure chart to GWLs of interest and relabel the specification
  df_linechart_exposure <- data_exposure %>%
    filter(warming_level %in% warming_levels_selected) %>%
    mutate(spec_label = factor(case_when(spec == "kotz" ~ "Impacts from all climate indicators\n(= Kotz et al., 2022)",
                                         spec == "status_quo" ~ "Impacts from ann. temp. only\n(= status quo)"),
                               levels = c("Impacts from ann. temp. only\n(= status quo)",
                                          "Impacts from all climate indicators\n(= Kotz et al., 2022)"))) %>%
    mutate(variable = if_else(variable == "imp_temp_Tonly_diff", "imp_total", variable)) %>%
    filter(variable == "imp_total")
  
  # make the linechart
  df_linechart_exposure %>%
    
    # NOTE: we flip the sign of perc05_threshold because impacts are labelled as 'damages'
    ggplot(aes(-perc05_threshold, pop_under_threshold)) +
    
    # one line for each warming level, different linetypes by specification
    geom_line(aes(linetype = spec_label, colour = paste0(warming_level, "\u00B0C"))) +

    # we add an errorbar for the +3C warming level & 15%GDP threshold, marking the difference between the specs
    annotate("errorbar", x = 0.15, width = 0.005, colour = "black", alpha = 0.9,
             ymin = df_linechart_exposure %>% filter(perc05_threshold == -0.15, warming_level == 3, spec == "status_quo") %>% pull(pop_under_threshold),
             ymax = df_linechart_exposure %>% filter(perc05_threshold == -0.15,  warming_level == 3, spec == "kotz") %>% pull(pop_under_threshold)) +
    
    # via annotated text label, we explain how to read the errorbar
    annotate("text", x = 0.155, y = 0.81, size = 7.5/.pt,
             label = paste0("Increase for risk of\n15%GDP loss by\nincluding all climate\nindicators: +",
                            100*round(df_linechart_exposure %>% filter(perc05_threshold == -0.15, warming_level == 3) %>%
                                        summarise(diff = pop_under_threshold[spec == "kotz"] - pop_under_threshold[spec == "status_quo"]) %>% pull(diff), 2), "%POP"),
             hjust = 0) +
    
    # set scales, labels and legend items
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::percent) +
    labs(linetype = NULL, y = "Global POP share of countries with\n>= 5% chance of damages above threshold",
         x = "Threshold (GDP loss)", colour = "Global warming level") +
    guides(colour = guide_legend(byrow = T),
           linetype = guide_legend(byrow = T)) +
    scale_linetype_manual(values = c("dotted", "solid")) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.spacing.y = unit(-0.2, 'cm'),
          legend.background = element_rect(fill = "transparent", colour = "transparent"))
}

# create the linechart for 1.5C, 2C and 3C
figure5a <- figure_exposure_linechart(df_pop_exposure, warming_levels_selected = c(1.5, 2, 3))

# write a function to create the population exposure table chart
figure_exposure_table <- function(data_exposure = NULL,
                                  warming_levels_selected = c(1.5, 2, 3)) {
  
  data_exposure %>%
    # rename imp_temp_Tonly_diff to immp_total
    mutate(variable = if_else(variable == "imp_temp_Tonly_diff", "imp_total", variable)) %>%
    
    # subset to threshold and GWLs to be shown in the table
    filter(variable %in% c("imp_total"),
           warming_level %in% warming_levels_selected,
           perc05_threshold %in% c(0, -0.05, -0.1, -0.15, -0.2)) %>%
    
    # clean up labels
    mutate(variable = rewrite_impact_label(variable),
           perc05_threshold_label = paste0(100*round(perc05_threshold, 2), "%"),
           spec_label = factor(case_when(spec == "kotz" ~ "Impacts from all climate indicators\n(= Kotz et al., 2022)",
                                         spec == "status_quo" ~ "Impacts from ann. temp. only\n (= status quo)"),
                               levels = c("Impacts from ann. temp. only\n (= status quo)",
                                          "Impacts from all climate indicators\n(= Kotz et al., 2022)")),
           warming_level_label = paste0(warming_level, "\u00B0C")) %>%
    
    # set up plot (NOTE: we order the character labels by the underlying numeric values)
    ggplot(aes(reorder(warming_level_label, warming_level), reorder(perc05_threshold_label, perc05_threshold))) +
    
    # add a tile coloured based on population falling under threshold
    geom_tile(aes(fill = pop_under_threshold)) +
    
    # one facet for main results and status quo results
    facet_wrap(~ spec_label) +
    
    # add the (rounded) population share under threshold to each tile
    geom_text(aes(label = paste0(100*round(pop_under_threshold, 2), "%")), size = 3) +
    
    # set color scales, labels and plot style
    scale_fill_gradient2(labels = scales::percent, high = "#F8766D", low = "white") +
    guides(fill = "none") +
    labs(x = "Global warming level", y = "Threshold (GDP loss)",
         subtitle = "Selected values of global population exposed to >=5% risk (from Panel a)") +
    theme_classic() +
    theme(plot.subtitle=element_text(size = 9, face = "bold"),
          axis.title = element_text(size = 9))
}

# make the table
figure5b <- figure_exposure_table(df_pop_exposure)

# create figure 5 and save out
figure5 <- ggarrange(figure5a,
          figure5b,
          heights = c(1.9, 1),
          labels = "auto",
          nrow = 2, ncol = 1)
ggsave(paste0("../sharepoint/Figures/",
              Sys.Date(), " Figure5.pdf"),
       plot = figure5,
       width = 5.2, height = 6)

# make an SI version for all six climate indicators using the Kotz et al specification
figure_exposure_table_siversion <- function(data_exposure = NULL, warming_levels_selected = c(1.5, 2, 3)) {
  
  # NOTE: steps are similar as the ones above
  
  data_exposure %>%
    filter(warming_level %in% warming_levels_selected) %>%
    filter(!variable %in% c("imp_nontemp", "imp_total", "imp_varextremes", "imp_temp", "imp_temp_Tonly_diff"),
           spec == "kotz",
           perc05_threshold > -0.05) %>%
    mutate(variable = rewrite_impact_label(variable)) %>%
    ggplot(aes(paste0(warming_level, "\u00B0C"), perc05_threshold )) + geom_tile(aes(fill = pop_under_threshold)) + facet_wrap(~ variable) +
    geom_text(aes(label = paste0(100*round(pop_under_threshold, 2), "%")), size = 3) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_gradient2(labels = scales::percent, high = "#F8766D", low = "white") +
    guides(fill = "none") +
    labs(x = "Global warming level", y = "Threshold (in %GDP)",
         subtitle = "% of global POP living in countries with >= 5% chance of GDP impacts worse than...") +
    theme_classic()
}

# create and save out
figure_exposure_table_siversion(df_pop_exposure)
ggsave(paste0("../sharepoint/Figures/",
              Sys.Date(), " SI_poptailrisk_perc05_bc_fulldistr_le.png"),
       width = 6.35, height = 4)

# clean up
rm(figure5a, figure5b, figure5, df_pop, df_pop_exposure, df_adm0_with_pop)
gc()


################################################################################
### Additional SI charts on variance decomposition  ############################
################################################################################

# plot the absolute variance using the Schwarzwald & Lenssen 2022 decomposition
df_var_bc_fulldistr_le_schwarzwald %>%
  figure_vardisaggr_global_schwarzwald(use_variance_share = F, data_fulldistr_aggr = df_global_bc_fulldistr_le_aggr,
                                       var_selection = c("imp_total", "imp_temp_diff", "imp_Pt_diff", "imp_varextremes",
                                                         "imp_Tstd_diff", "imp_Wn_diff", "imp_wet_days_1_diff", "imp_vwet_days1_am_99p9_diff"))
ggsave(file.path("..", "sharepoint", "Figures", paste0(Sys.Date(), " FigurSI_var_schwarzwald_absolute.png")),
       height = 4, width = 6.5)

# plot the coefficient of variation by GWL & climate indicator
df_global_bc_fulldistr_le_aggr %>%
  # select GWL as well as means and variances for each impact channel
  select(warming_level, ends_with("_var"), ends_with("_mean")) %>%
  # discard status quo impacts
  select(-contains("Tonly")) %>%
  # convert variances to st.dv.
  mutate_at(vars(ends_with("_var")), sqrt) %>%
  # convert to LONG format, create separate column for summary stat (mean or st.dv.), clean up the variable name
  pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "value") %>%
  mutate(summarystat = case_when(str_detect(variable, "_var$") ~ "stdv",
                                 str_detect(variable, "_mean$") ~ "mean")) %>%
  mutate(variable = variable %>% str_remove("_var$") %>% str_remove("_mean$")) %>%
  # convert to WIDE (separate columns for mean and st.dv.), then calculate coefficient of variation and clean up climate indicator labels
  pivot_wider(names_from = "summarystat", values_from = "value") %>%
  mutate(coef_variation = stdv/abs(mean),
         variable = rewrite_impact_label_linebreaks(variable)) %>%
  # subset to GWLs of interest
  filter(warming_level %in% c(1.5, 2, 3)) %>%
  # make bar chart
  ggplot(aes(paste0(warming_level, "\u00B0C"), coef_variation)) + geom_col(alpha = 0.7, fill = "#619CFF") +
  # one facet per climate indicator
  facet_wrap(~ variable, nrow = 2, scales = "free_y") +
  labs(y = "Coefficient of variation\nof global GDP impacts", x = "Global warming level") +
  theme_classic()
ggsave(file.path("..", "sharepoint", "Figures", paste0(Sys.Date(), " FigureSI_coeficient_variation.png")),
       height = 3, width = 6)


### Variance decomposition following Hsiang et al 2017

# plot variances
figure_vardisaggr_global_hsiang <- function(data = NULL,
                                            warming_levels_shown = c(1.5, 2, 3),
                                            nrow_used = NA, ncol_used = NA,
                                            show_interaction = TRUE, interaction_as_bar = TRUE,
                                            absolute_value_rescaling = FALSE,
                                            model_label = NULL) {
  data %>%
    filter(warming_level %in% warming_levels_shown) %>%
    # discard the temperature damages using the 'status quo' approach
    filter(!variable == "imp_temp_Tonly_diff") %>%
    # if absolute_value_rescaling is TRUE, we rescale the variance if interaction term is negative and then take absolute values
    {if(absolute_value_rescaling) mutate(.,
                                         variance_total = if_else(interaction < 0, variance_total - 2*interaction, variance_total),
                                         interaction = abs(interaction)) else .} %>%
    {if(interaction_as_bar) pivot_longer(., cols = c("monte_carlo_draw", "model", "scenarioensembleyear", "interaction"),
                                         names_to = "variance_source", values_to = "variance") else pivot_longer(., cols = c("monte_carlo_draw", "model", "scenarioensembleyear"),
                                                                                                                 names_to = "variance_source", values_to = "variance")
    } %>%
    group_by(GID_0, warming_level, variable) %>%
    # if we show the interaction term, we divide variance components by total variance - otherwise by the sum of variance components ('variance_explained')
    {if(interaction_as_bar) mutate(., variance_share = variance/variance_total) else mutate(., variance_share = variance/variance_explained)} %>%
    mutate(variance_total_relative = variance_total/variance_explained) %>%
    # discard the interaction term if we do not show it
    {if(!interaction_as_bar) filter(., variance_source != "interaction") else . } %>%
    # clean up the variance component labels
    mutate(variable = rewrite_impact_label_linebreaks(variable),
           variance_source = variance_source %>% str_replace("^model$", "Climate model\nuncertainty") %>%
             str_replace("^year$", "Internal\nvariability") %>%
             str_replace("^scenarioensembleyear$", "Internal\nvariability") %>%
             str_replace("^monte_carlo_draw$", "Dose-response function\nuncertainty") %>%
             str_replace("^scenario$", "Scenario\nuncertainty") %>%
             str_replace("^interaction$", "Interaction/\nshared variance")) %>%
    # convert to factor to fix the item order
    mutate(variance_source = factor(variance_source, levels = c("Dose-response function\nuncertainty", "Internal\nvariability", "Climate model\nuncertainty",
                                                                "Scenario\nuncertainty", "Interaction/\nshared variance"))) %>%
    ungroup() %>%
    
    # make bar chart
    ggplot(aes(paste0(warming_level, "\u00B0C"), variance_share)) +
    
    # if we show the interaction, make a horizontal line at 100% (since sum of bars can exceed 100%)
    {if(show_interaction) geom_hline(yintercept = 1, alpha = 0.3, linetype = "dashed")} +
    
    # stacked bars for variance components
    geom_bar(aes(fill = variance_source), stat = "identity", position = position_stack(vjust = 0.5, reverse = F), alpha = 0.7) +
    
    # label variance shares
    geom_label(aes(fill = variance_source,
                   label = paste0(100*(round(variance_share, ifelse(round(variance_share, 2) == 0, 3, 2))), "%")),
               position = position_stack(vjust = 0.5, reverse = F),
               size = 2.5,
               show.legend = F) +
    
    # if we do not show interactions as separate bars, add a point for the total variance
    {if(show_interaction & !interaction_as_bar) geom_point(aes(y = variance_total_relative, shape = "Total variance\nplus interaction"), colour = "black", alpha = 0.5)} +
    
    # set labels, scales and plot style
    labs(x = "Global warming level", y = "Share in variance", colour = NULL, fill = NULL, shape = NULL) +
    {if(!is.null(model_label)) labs(subtitle = paste0("Central model: ", model_label)) } +
    scale_fill_manual(values = c("#F8766D", "#619CFF", "#00BA38", "grey")) +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(~ variable, nrow = 2, ncol = 4, scales = ifelse(show_interaction & !absolute_value_rescaling, "free_y", "fixed")) +
    theme_classic() +
    theme(legend.position = "bottom")
}

# write a function to set up the tibble for variance decomposition charts a la Hsiang et al 2017 for different median-like models
load_vardisaggr_hsiang <- function(model_selected = NULL) {
  
  # throw an error if marginal variances using the model selected have not been calculated yet
  if(!file.exists(file.path("data", paste0("df_global_bc_fulldistr_le_hsiangvar_mc_", model_selected, ".feather")))) stop(paste0("Marginal variance of dose-response draws missing for model ", model_selected))
  if(!file.exists(file.path("data", paste0("df_global_bc_fulldistr_le_hsiangvar_internal_", model_selected, ".feather")))) stop(paste0("Marginal variance for internal variability missing for model ", model_selected))
  
  # otherwise, load the respective marginal variances, merge in total variance and calculate the residual variance
  df_out <- bind_rows(
    # load the marginal variance for model uncertainty - this is independent of the median-like model chosen
    file.path("data", "df_global_bc_fulldistr_le_hsiangvar_model.feather") %>% read_feather(),
    # load MC dose-response & internal variability variances - these depend on which CMIP6 model was chosen as central value
    file.path("data", paste0("df_global_bc_fulldistr_le_hsiangvar_mc_", model_selected, ".feather")) %>% read_feather(),
    file.path("data", paste0("df_global_bc_fulldistr_le_hsiangvar_internal_", model_selected, ".feather")) %>% read_feather()
  ) %>% pivot_wider(values_from = "variance", names_from = "variance_source", -n) %>%
    # merge in the true total variance from the summary stats tibble, so we can calculate the residual 
    left_join(df_global_bc_fulldistr_le_aggr %>% select(warming_level, ends_with("_var")) %>%
                pivot_longer(cols = ends_with("_var"), names_to = "variable", values_to = "variance_total") %>%
                mutate(GID_0 = "Global") %>%
                mutate(variable = str_remove(variable, "_var$")),
              by = c("GID_0", "warming_level", "variable")) %>%
    # calculate the residual
    mutate(variance_explained = monte_carlo_draw + model + scenarioensembleyear,
           interaction = variance_total - variance_explained)
  
  return(df_out)
}

# plot the average GDP impact by model to identify (near-)median models
ggarrange(
  plotlist = map(c("imp_total", "imp_temp_diff", "imp_Pt_diff", "imp_varextremes",
                   "imp_Tstd_diff", "imp_Wn_diff", "imp_wet_days_1_diff", "imp_vwet_days1_am_99p9_diff"),
                 .f  = function(varname) {
                   df_global_bc_fulldistr_le_gwl %>%
                     # subset to +3C GWL, select impact variable of interest and average
                     filter(warming_level == 3) %>% group_by(model) %>%
                     select(warming_level, model, !!as.name(varname)) %>%
                     summarise_at(vars(starts_with("imp")), mean) %>%
                     # convert to LONG format
                     pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "imp") %>%
                     # clean up climate indicator names
                     mutate(variable = rewrite_impact_label(variable)) %>%
                     # for each climate indicator, order by mean impact per model and identify the median model
                     group_by(variable) %>%
                     arrange(variable, imp) %>% 
                     mutate(ind = 1:n(),
                            n = n()) %>%
                     mutate(is_median_model = if_else(n %% 2 == 0, ind == n/2, ind == n/2 + 0.5)) %>%
                     ungroup() %>%
                     
                     # make a bar chart, order models by their mean impact, and color the median model differently
                     ggplot(aes(reorder(model, ind), imp)) + geom_col(aes(fill = is_median_model), alpha = 0.7) +
                     scale_y_continuous(labels = scales::percent, n.breaks = 3) +
                     scale_fill_manual(values = c("#619CFF", "#F8766D")) +
                     labs(x = NULL, y = paste0("Mean global GDP\nimpact at +3\u00B0C"), fill = NULL) +
                     guides(fill = "none") +
                     coord_flip() +
                     facet_wrap(~ variable, scales = "free")
                 }), nrow = 2, ncol = 4)
ggsave(file.path("..", "sharepoint", "Figures", paste0(Sys.Date(), " FigureSI_medianmodel_GWL3.png")),
       height = 8, width = 12)

# create the chart for different central models and save it out
for(model_jj in c("KACE-1-0-G", "EC-Earth3-Veg", "CESM2", "MPI-ESM1-2-LR")) {
  
  load_vardisaggr_hsiang(model_jj) %>% figure_vardisaggr_global_hsiang(nrow_used = 2, ncol_used = 4,
                                                                           model_label = model_jj)
  ggsave(file.path("..", "sharepoint", "Figures", paste0(Sys.Date(), " Figure2_var_hsiang_", model_jj, ".png")),
         height = 7, width = 10)
  
}



################################################################################
### SI CHARTS ON INCLUDING Tx5d BASED ON Callahan & Mankin 2022 ################
################################################################################

# import summary statistics for including Tx5d (based on point estimates for dose-rseponse function parameters)
df_global_bc_pointdistr_le_tx5d_aggr <- read_feather(file.path("data", "df_global_bc_pointdistr_le_tx5d_aggr.feather"))

# compare results incl. Tx5d to our overall results (using large ensemble & full Monte Carlos)
bind_rows(df_global_bc_pointdistr_le_tx5d_aggr %>% mutate(label = "Including Tx5d"),
          df_global_bc_fulldistr_le_aggr %>% mutate(label = "Excluding Tx5d\n(= main results)")) %>%
  select(warming_level, ends_with("mean"), label) %>%
  filter(warming_level %in% c(3)) %>%
  pivot_longer(cols = ends_with("mean"), names_to = "variable", values_to = "imp") %>%
  # discard aggregates and status quo impacts, which are not really meaningful here
  filter(!variable %in% c("imp_total_mean", "imp_varextremes_mean")) %>%
  mutate(alpha_category = case_when(variable != "imp_temp_Tonly_diff_mean" ~ "1",
                                    str_detect(label, "Excluding") ~ "2",
                                    TRUE ~ "3"),
         variable = variable %>% str_remove("_mean") %>% rewrite_impact_label(),
         imp = replace_na(imp, 0)) %>%
  
  # make a bar chart for mean global GDP impacts by impact channel
  ggplot(aes(variable, imp)) +
  geom_col(aes(fill = label, alpha = alpha_category), position = position_dodge()) +
  
  # flip axes and set scales, labels, plot style
  coord_flip() + 
  scale_alpha_manual(values = c(1, 0.2, 0)) +
  labs(y = "Global GDP impact at +3\u00B0C", x = NULL, fill = NULL) +
  guides(alpha = "none") +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom")

ggsave(file.path("..", "sharepoint", "Figures", 
                 paste0(Sys.Date(), " SI_tx5d_globalresults_pointestimates_gwl3.png")),
       width = 4.25, height = 3.5)


################################################################################
##################### SI CHARTS ON LARGE ENSEMBLE INCLUSION ####################
################################################################################

# NOTE: as done above, we convert values to basis points for readability
rescale_imp_factor <- 10^4



### Part a) assess internal variability of large ensembles


### extract large ensemble values from distribution without dose-response function variation ('pointdistr') and calculate variances (= internal variability)

# CESM2-LE
df_var_cesm_pointdistr <- df_global_bc_pointdistr_le_gwl %>%
  filter(model == "CESM2", scenario == "ssp370") %>%
  group_by(model, scenario, warming_level) %>%
  # rescale impacts and calculate variance
  mutate_at(vars(starts_with("imp")), ~ .x * rescale_imp_factor) %>%
  summarise_at(vars(starts_with("imp")), var)

# MPI-ESM1-2-LR
df_var_mpi_pointdistr <- df_global_bc_pointdistr_le_gwl %>%
  filter(model == "MPI-ESM1-2-LR", scenario == "ssp370") %>%
  group_by(model, scenario, warming_level) %>%
  # rescale impacts and calculate variance
  mutate_at(vars(starts_with("imp")), ~ .x * rescale_imp_factor) %>%
  summarise_at(vars(starts_with("imp")), var)

# calculate variance for each model run ('ensemble') = inter-annual variability
df_var_cesm_pointdistr_by_ensemble <- df_global_bc_pointdistr_le_gwl %>%
  filter(model == "CESM2", scenario == "ssp370") %>%
  group_by(model, scenario, ensemble, warming_level) %>%
  mutate_at(vars(starts_with("imp")), ~ .x * rescale_imp_factor) %>%
  summarise_at(vars(starts_with("imp")), var)

# repeat for MPI-ESM1-2-LR
df_var_mpi_pointdistr_by_ensemble <- df_global_bc_pointdistr_le_gwl %>%
  filter(model == "MPI-ESM1-2-LR", scenario == "ssp370") %>%
  group_by(model, scenario, ensemble, warming_level) %>%
  mutate_at(vars(starts_with("imp")), ~ .x * rescale_imp_factor) %>%
  summarise_at(vars(starts_with("imp")), var)

# function to make comparison chart for internal variability
compare_internal_variability <- function(data_by_ensemble = NULL, # internal variability by ensemble member
                                         data_full_ensemble = NULL # internal variability of the full ensemble
                                         ) {
  
  data_by_ensemble %>%
    # subset to GWLs of interest
    filter(warming_level %in% c(1.5, 2, 3)) %>%
    
    # grab impacts of interest and convert to LONG format
    pivot_longer(cols = c(imp_total, imp_temp_diff, imp_varextremes),
                 names_to = "variable", values_to = "internal_variability") %>%
    
    # clean up impact channel labels
    mutate(variable = rewrite_impact_label(variable)) %>%
    
    # initiate plot
    ggplot(aes(paste0("+", warming_level, "C"), internal_variability)) +
    
    # show internal variability of each ensemble member as jitter plot
    geom_jitter(aes(y = internal_variability), colour = "grey", alpha = 0.3, width = 0.2) +
    
    # boxplot without whiskers for lower/upper quartiles & median
    geom_boxplot(aes(group = warming_level), outlier.shape = NA, coef = 0, fill = NA) +
    
    # add average across ensemble members as a point
    stat_summary(
      geom = "point",
      fun = "mean",
      aes(col = "Average across ensemble members",
          shape = "Average across ensemble members"),
      size = 3, alpha = 0.8
    ) +
    
    # bring the internal variability of the full ensemble into the same format and show as a different point
    geom_point(data = data_full_ensemble %>% filter(warming_level %in% c(1.5, 2, 3)) %>%
                 pivot_longer(cols = c(imp_total, imp_temp_diff, imp_varextremes),
                              names_to = "variable", values_to = "internal_variability") %>%
                 mutate(variable = rewrite_impact_label(variable),
                        is_large_ensemble = T),
               aes(y = internal_variability,
                   shape = "Entire large ensemble",
                   colour = "Entire large ensemble"),
               size = 3, alpha = 0.8) +
    
    # one facet per impact channel
    facet_wrap(~ variable, scales = "free_y") +
    
    # set labels
    labs(colour = NULL, shape = NULL, x = "Global warming level", y = "Internal variability of\n global GDP impacts\n(in squared bps)") +
    theme(legend.position = "bottom")
}

# make the 2x1 grid chart (one row for CESM2-LE and the MPI ensemble)
ggarrange(compare_internal_variability(df_var_cesm_pointdistr_by_ensemble, df_var_cesm_pointdistr) +
            labs(title = "CESM2-LE", subtitle = "One data point per ensemble member"),
          compare_internal_variability(df_var_mpi_pointdistr_by_ensemble, df_var_mpi_pointdistr) +
            labs(title = "MPI-ESM1-2-LR large ensemble (30 realizations only)",
                 subtitle = "One data point per ensemble member"),
          nrow = 2)

# save out
ggsave(file.path("..", "sharepoint", "Figures", 
                 paste0(Sys.Date(), " SI_largeensemble_internalvar.png")),
       width = 6, height = 6)

# do a variance decomposition by treating ensemble members as models for MPI-ESM1-2-LR
df_varcomp_ensemble_mpi <- df_global_bc_pointdistr_le_gwl %>%
  
  # subset to the large ensemble for SSP3-7.
  filter(model == "MPI-ESM1-2-LR", scenario == "ssp370") %>%
  
  # grab variables of interest and discard the status quo impacts ('imp_temp_Tonly_diff')
  select(model, scenario, ensemble, year, warming_level, imp_total, imp_varextremes, ends_with("diff"),
         -contains("Tonly")) %>%
  
  # set to GWLs of interest
  filter(warming_level %in% c(1.5, 2, 3)) %>%
  
  # convert to LONG format and rescale impacts
  pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "imp") %>%
  mutate(imp = imp * rescale_imp_factor) %>%
  
  # for each model run, warming_level and impact channel, calculate conditional mean and variance
  group_by(ensemble, warming_level, variable) %>%
  summarise(var_within_ensemble = var(imp),
            mean_within_ensemble = mean(imp)) %>% ungroup() %>%
  
  # across model runs, calculate variance of conditional mean (= between-run variance) and mean of conditional variance (= within-run)
  group_by(warming_level, variable) %>%
  summarise(ensemble_uncertainty = var(mean_within_ensemble),
            interannual = mean(var_within_ensemble))

# make a variance decomposition chart
df_varcomp_ensemble_mpi %>%
  
  # convert to long format, so we have one row each for ensemble uncertainty and inter-annual variability
  pivot_longer(cols = c(ensemble_uncertainty, interannual)) %>%
  
  # calculate shares
  group_by(warming_level, variable) %>% mutate(value = value/sum(value)) %>% ungroup() %>%
  
  # clean up labels
  mutate(name = case_when(name == "ensemble_uncertainty" ~ "Variance between ensemble runs",
                          name == "interannual" ~ "Variance within ensemble runs")) %>%
  mutate(variable = rewrite_impact_label_linebreaks(variable)) %>%
  
  # make the bar chart
  ggplot(aes(paste0(warming_level, "\u00B0C"), value)) +
  geom_col(aes(fill = name)) +
  
  # one facet per impact channel
  facet_wrap(~ variable, nrow = 2) +
  
  # set labels and scales
  labs(fill = NULL, x = "Global warming level", y = "Share in variance") +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom")

ggsave(file.path("..", "sharepoint", "Figures", 
                 paste0(Sys.Date(), " SI_largeensemble_vardecomposition_MPI-ESM1-2-LR.png")),
       width = 6, height = 4)


# clean up the environment
rm(df_var_cesm_pointdistr, df_var_cesm_pointdistr_by_ensemble,
   df_var_mpi_pointdistr, df_var_mpi_pointdistr_by_ensemble, df_varcomp_ensemble_mpi)
gc()



### Part b) assess tails of large ensembles vs single-run distribution of many models

# NOTE: summarystats_for_distribution() requires weights but for distributions involving
# only the same large ensemble (either single run or all runs), this will be uniform
# and hence it does not matter which weights are provided

# calculate summary stats for CESM2-LE
df_global_bc_pointdistr_le_aggr_cesm2 <- df_global_bc_pointdistr_le_gwl %>%
  filter(model == "CESM2", scenario == "ssp370") %>%
  select(model, scenario, ensemble, year, warming_level,
         imp_total, imp_varextremes, ends_with("diff")) %>%
  summarystats_for_distribution(gcm_weights_used = df_gcm_weights_le)

# summary stats for only first run of CESM2-LE
df_global_bc_pointdistr_le_aggr_cesm2_singlerun <- df_global_bc_pointdistr_le_gwl %>%
  filter(model == "CESM2", scenario == "ssp370", ensemble == "r1i1p1f1") %>%
  select(model, scenario, ensemble, year, warming_level,
         imp_total, imp_varextremes, ends_with("diff")) %>%
  summarystats_for_distribution(gcm_weights_used = df_gcm_weights_le)

# summary stats of MPI-ESM1-2-LR ensemble
df_global_bc_pointdistr_le_aggr_mpi <- df_global_bc_pointdistr_le_gwl %>%
  filter(model == "MPI-ESM1-2-LR", scenario == "ssp370") %>%
  select(model, scenario, ensemble, year, warming_level,
         imp_total, imp_varextremes, ends_with("diff")) %>%
  summarystats_for_distribution(gcm_weights_used = df_gcm_weights_le)

# summary stats for only first run of MPI ensemble
df_global_bc_pointdistr_le_aggr_mpi_singlerun <- df_global_bc_pointdistr_le_gwl %>%
  filter(model == "MPI-ESM1-2-LR", scenario == "ssp370", ensemble == "r1i1p1f1") %>%
  select(model, scenario, ensemble, year, warming_level,
         imp_total, imp_varextremes, ends_with("diff")) %>%
  summarystats_for_distribution(gcm_weights_used = df_gcm_weights_le)

# summarystats for single-run distribution where we subset to SSP3-7.0 and discard all but
# the first runs of the two large ensembles
# NOTE: since this distribution features different models, we need weights. However,
# as we show in the next step, weights are uniform for +3C (because each model
# reaches this GWL only for one scenario, SSP3-7.0, which we use here)
df_comp_aggr <- df_global_bc_pointdistr_le_gwl %>%
  # subset to SSP3-7.0
  filter(scenario == "ssp370") %>%
  # discard the large ensembles
  filter(!(model == "CESM2" & ensemble != "r1i1p1f1")) %>%
  filter(!(model == "MPI-ESM1-2-LR" & ensemble != "r1i1p1f1")) %>%
  select(model, scenario, ensemble, year, warming_level,
         imp_total, imp_varextremes, ends_with("diff")) %>%
  summarystats_for_distribution(gcm_weights_used = df_gcm_weights_cmip6)

# show that model weights are uniform
df_gcm_weights_cmip6 %>% filter(warming_level == 3) %>% summary()
# -> means that in the chart below, we can take simple average and result is valid

# figure to compare the distributions
figure_distribution_largeensemble <- function(data = NULL, # data for the full distribution
                                              data_aggr = NULL, # summary stats for distribution across many models
                                              data_single_run_aggr = NULL, # summary stats for single run of ensemble
                                              data_fullensemble_aggr = NULL, # summary stats for entire ensemble
                                              single_run_id = "r1i1p1f1", # model run used for single run of ensembles
                                              model_label = NULL, # model label
                                              allmodel_label = "All models\n(Single runs)", # label for the many-models distribution
                                              imp_selected = "imp_total", # impact channel displayed
                                              warming_level_selected = c(3) # GWL displayed
                                              ) {
  
  # create labels for the different distribution that include the user-specific model_label
  single_run_label <- paste0(model_label, "\nSingle run")
  full_ensemble_label <- paste0(model_label, ifelse(str_detect(model_label, "MPI"),
                                                    "\nPartial ensemble", "\nFull ensemble"))
  
  # bind the relevant part of the distribution together and label them
  bind_rows(# a) the distribution for the single run of the ensemble
            data %>% filter(model == model_label, scenario == "ssp370", ensemble == single_run_id) %>% mutate(ensemble_label = single_run_label),
            # b) the distribution for the entire large ensemble
            data %>% filter(scenario == "ssp370", model == model_label) %>% mutate(ensemble_label = full_ensemble_label),
            # c) the many-model distribution with only first runs of each ensemble
            data %>% filter(scenario == "ssp370") %>%
              # discard the large ensembles
              filter(!(model == "CESM2" & ensemble != "r1i1p1f1")) %>%
              filter(!(model == "MPI-ESM1-2-LR" & ensemble != "r1i1p1f1")) %>%
              mutate(ensemble_label = allmodel_label)
  ) %>%
    
    # subset to GWLs of interest 
    filter(warming_level %in% warming_level_selected) %>%
    mutate(warming_level = paste0(warming_level, "\u00B0C")) %>%
    group_by(ensemble_label) %>%
    
    # rename the variable selected to y_var
    rename(y_var = imp_selected) %>%
    ggplot(aes(x = ensemble_label)) +
    
    # jitter plot for the entire distribution, by ensemble_label
    geom_point(aes(y = y_var, colour = ensemble_label),
               position = position_jitter(w = .15),
               size = 0.5,
               alpha = 0.1) +
    
    # grey triangle for distribution mean
    stat_summary(aes(y = y_var, colour = ensemble_label),
                 geom = "point", fun.y = "mean", size = 4, col = "darkgrey",
                 shape = 17) +
    
    # boxplot for summary stats by distribution
    geom_boxplot(data = bind_rows(data_aggr %>% mutate(ensemble_label = allmodel_label),
                                  data_single_run_aggr %>% mutate(ensemble_label = single_run_label),
                                  data_fullensemble_aggr %>% mutate(ensemble_label = full_ensemble_label)) %>%
                   # subset to variables of interest and convert to LONG format
                   select(warming_level, contains(imp_selected), ensemble_label) %>%
                   filter(warming_level %in% warming_level_selected) %>%
                   mutate(warming_level = paste0(warming_level, "\u00B0C")) %>%
                   pivot_longer(cols = starts_with("imp"), names_to = "summarystat", values_to = "imp") %>%
                   # clean up the summarystat column
                   mutate(summarystat = str_remove(summarystat, imp_selected) %>%  str_remove("^_")) %>%
                   group_by(ensemble_label, warming_level) %>%
                   # convert to WIDE so we have one column per summary stat (one row per climate indicator and distribution)
                   pivot_wider(names_from = "summarystat", values_from = "imp"),
                 
                 # choose boxplot statistics manually
                 stat = "identity",
                 aes(lower = perc25, upper = perc75,
                     middle = perc50, ymin = perc05, ymax = perc95, group = ensemble_label,
                     fill = ensemble_label),
                 width = .25,
                 # no dots for outliers
                 outlier.shape = NA,
                 alpha = 0.1) +
    
    # add a one-sided violin chart for the distributions
    geom_flat_violin(aes(y = y_var, fill = ensemble_label),
                     position = position_nudge(x = .2),
                     alpha = 0.25,
                     adjust = 0.5) +
    
    # one facet per impact channel
    facet_wrap(~ rewrite_impact_label(imp_selected)) +
    
    # set scales and labels
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c( "#619CFF", "#00BA38", "#F8766D")) +
    scale_color_manual(values = c( "#619CFF", "#00BA38", "#F8766D")) +
    scale_x_discrete(limits = rev) +
    labs(x = NULL,
         y = paste0("Global GDP impact at +", warming_level_selected, "\u00B0C")) +
    theme_classic() +
    theme(legend.position = "none")
}

# make the same chart for eight impact channels of interest and combine in a 4x2 grid chart
ggarrange(plotlist = map(c("imp_total", "imp_temp_diff", "imp_Pt_diff", "imp_varextremes",
                       "imp_Tstd_diff", "imp_Wn_diff", "imp_wet_days_1_diff", "imp_vwet_days1_am_99p9_diff"),
                        ~ figure_distribution_largeensemble(data = df_global_bc_pointdistr_le_gwl,
                                                            data_aggr = df_comp_aggr,
                                                            data_single_run_aggr = df_global_bc_pointdistr_le_aggr_cesm2_singlerun,
                                                            data_fullensemble_aggr = df_global_bc_pointdistr_le_aggr_cesm2,
                                                            model_label = "CESM2",
                                                            imp_selected = .x)
), nrow = 2, ncol = 4)

ggsave(file.path("..", "sharepoint", "Figures", 
                 paste0(Sys.Date(), " SI_largeensemble_tails_CESM2.png")),
       width = 13, height = 5.5)

# make the same chart for MPI
ggarrange(
  plotlist = map(c("imp_total", "imp_temp_diff", "imp_Pt_diff", "imp_varextremes",
                   "imp_Tstd_diff", "imp_Wn_diff", "imp_wet_days_1_diff", "imp_vwet_days1_am_99p9_diff"),
                 ~ figure_distribution_largeensemble(data = df_global_bc_pointdistr_le_gwl,
                                                     data_aggr = df_comp_aggr,
                                                     data_single_run_aggr = df_global_bc_pointdistr_le_aggr_mpi_singlerun,
                                                     data_fullensemble_aggr = df_global_bc_pointdistr_le_aggr_mpi,
                                                     model_label = "MPI-ESM1-2-LR",
                                                     imp_selected = .x)
  ), nrow = 2, ncol = 4)

ggsave(file.path("..", "sharepoint", "Figures", 
                 paste0(Sys.Date(), " SI_largeensemble_tails_MPI-ESM1-2-LR.png")),
       width = 15, height = 5.5)

# clean up the environment
rm(df_comp_aggr, df_global_bc_pointdistr_le_aggr_mpi_singlerun, df_global_bc_pointdistr_le_aggr_mpi,
   df_global_bc_pointdistr_le_aggr_cesm2_singlerun, df_global_bc_pointdistr_le_aggr_cesm2)
gc()



################################################################################
########################## ALTERNATIVE BIAS CORRECTION #########################
################################################################################

# calculate summary stats for the single model run using alternative bias correction
df_global_aggr_altbiascor <- file.path(dir_imp, "pointestimates_stagg_altbiascor", "global_.feather") %>%
  read_feather() %>%
  left_join(df_gwl, by = c("model", "scenario", "ensemble", "year")) %>%
  filter(warming_level %in% c(1.5, 2, 3)) %>%
  select(model, scenario, ensemble, year, warming_level, imp_total, imp_varextremes, ends_with("diff")) %>%
  summarystats_for_distribution(gcm_weights_used = df_gcm_weights_le,
                                is_input_datatable = F,
                                is_input_adm0 = F,
                                warming_levels_used = c(1.5, 2, 3))

# calculate summary stats just for the corresponding MPI-LR single run (r1i1p1f1, SSP3-7.0)
df_global_mpi_run_aggr <- df_global_bc_fulldistr_le_gwl %>%
  filter(model == "MPI-ESM1-2-LR" & scenario == "ssp370" & ensemble == "r1i1p1f1") %>%
  select(model, scenario, ensemble, year, warming_level, imp_total, imp_varextremes, ends_with("diff")) %>%
  summarystats_for_distribution(gcm_weights_used = df_gcm_weights_le,
                                is_input_datatable = F,
                                is_input_adm0 = F,
                                warming_levels_used = c(1.5, 2, 3))

# compare global results at +3C
bind_rows(df_global_mpi_run_aggr %>% mutate(label = "Main results\n(= bias-correcting\nannual indicators)"),
          df_global_aggr_altbiascor  %>% mutate(label = "Bias-correcting\ndaily temp. & precip.")) %>%
  
  # subset to variables of interest
  select(warming_level, ends_with("mean"),
         ends_with("perc10"), ends_with("perc90"), label) %>%
  
  # convert to LONG format and create separate columns for summarystat & impact channel
  pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "imp") %>%
  mutate(summarystat = str_extract(variable, "(?<=_)[:alnum:]+$"),
         variable = str_remove(variable, "_[:alnum:]+$")) %>%
  
  # convert to WIDE format such that we have one column per summarystat and one row per impact channel
  pivot_wider(names_from = "summarystat", values_from = "imp") %>%
  
  # discard aggregate for variability and extremes
  filter(!variable %in% c("imp_varextremes")) %>%
  mutate(variable = rewrite_impact_label(variable),
          alpha_label = str_detect(variable, "status quo")) %>%
  
  # make a bar chart with mean GDP impacts
  ggplot(aes(rewrite_impact_label(variable))) +
  geom_col(aes(y = mean, fill = label, alpha = alpha_label), position = position_dodge2()) +
  
  # one facet per GWL
  facet_wrap(~ paste0("+", warming_level, "\u00B0C"), scales = "free_x") +
  
  # set labels, scales and legend items
  labs(y = "Global GDP impact at +3\u00B0C", x = NULL, fill = NULL, alpha = NULL,
       subtitle = "Example model run (MPI-ESM1-2-LR, SSP3-7.0, r1i1p1f1)") +
  guides(alpha = "none") +
  scale_y_continuous(labels = scales::percent) +
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(legend.position = "bottom") +
  coord_flip()

ggsave(file.path("..", "sharepoint", "Figures", 
                 paste0(Sys.Date(), " SI_altbiascor_MPI-ESM1-2-LR.png")),
       width = 6.5, height = 3.5)

# clean up the environment
rm(df_global_mpi_run_aggr, df_global_aggr_altbiascor)
gc()



################################################################################
######################### SUMMARY STATS IN THE TEXT ############################
################################################################################

### Abstract

# number of models
unique(df_global_bc_fulldistr_le_gwl$model) %>% length()

# worst ADM0-level mean impact at +3C
df_adm0_bc_fulldistr_le_aggr %>% filter(warming_level == 3, variable == "imp_total") %>%
  slice_min(mean, n = 5) %>%
  select(warming_level, GID_0, mean) %>% mutate(mean = round(mean, 2))


### Main: Projecting GDP impacts for precipitation and temperature indicators

# show that +3C is only reached for SSP3-7.0 across (as implied by Figure 1) 
df_gwl %>% filter(warming_level == 3) %>% count(scenario)


### Main: Global Results

# NOTE: global summary stats by indicators (mean, upper/lower deciles) can be found in the respective SI table

# combined impact of variability & extremes at +3C (not covered in the table)
df_global_bc_fulldistr_le_aggr %>% filter(warming_level == 3) %>%
  select(warming_level, contains("varextr"))

# agreement of the contribution on mean sign
df_global_bc_fulldistr_le_aggr %>% filter(warming_level == 3) %>%
  select(warming_level, contains("agree")) %>%
  pivot_longer(cols = starts_with("imp"))



### Main: Country-level Results

# who are the countries with the highest certainty for extreme precip damages in Figure 3?
df_adm0_bc_fulldistr_le_aggr %>% filter(warming_level == 3, variable == "imp_vwet_days1_am_99p9_diff") %>%
  slice_max(n=5, order_by = agreement_meansign)

# show that no sovereign country in the sample has positive total impact
df_adm0_bc_fulldistr_le_aggr %>% filter(warming_level == 3) %>%
  filter(variable == "imp_total") %>% 
  slice_max(mean, n = 5)

# show that annual precipitation benefits most countries on average
df_adm0_bc_fulldistr_le_aggr %>% filter(warming_level == 3, variable == "imp_Pt_diff") %>%
  summarise(sum(mean > 0))


### Main: Overall impact of including variability and extremes

# calculate the difference in marginal effects for every initial temperature between -10 and +30C
# NOTE: this uses the marginal_effect_temperature() defined above for Figure 3
summary(100*(marginal_effect_temperature(T_0 = seq(-10, 30, by = 0.01), coef_used = get_coefficients()[, 1]) -
               marginal_effect_temperature(T_0 = seq(-10, 30, by = 0.01), coef_used = get_coefficients(model_used = "T_Pt")[, 1])))

# which countries have lower damages for status quo
df_adm0_bc_fulldistr_le_aggr  %>% filter(variable %in% c("imp_total", "imp_temp_Tonly_diff")) %>%
  select(warming_level, GID_0, variable, mean) %>%
  pivot_wider(names_from = "variable", values_from = "mean") %>%
  mutate(diff_statusquo = imp_total - imp_temp_Tonly_diff) %>%
  arrange(desc(diff_statusquo))


### Methods

# number of model-scenario combinations (counting one realization per model)
df_global_bc_fulldistr_le_gwl %>% count(model, scenario) %>% nrow()
# number of model-scenario-realization combinations
df_global_bc_fulldistr_le_gwl %>% count(model, scenario, ensemble) %>% nrow()

# sovereign countries not covered in the SSP database
df_ssp <- readRDS("data/df_ssp.rds")
sovereign_countries_gid0[!sovereign_countries_gid0 %in% df_ssp$GID_0]
