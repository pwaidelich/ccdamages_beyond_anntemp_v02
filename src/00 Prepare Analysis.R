# clean the environment
rm(list = ls())

# load packages
library(tidyverse)  # for general data wrangling
library(readxl)     # for importing Excel workbook data
library(xtable)     # for exporting Latex tables
library(rnaturalearth)   # for geodata to create map charts
library(rnaturalearthdata) # for low-resolution ADM0-level shapefiles

# map the directory for the ADM1-level climate data
dir_climate_data_adm1 <- file.path("/net", "cfc", "landclim1", "pwaidelich", "cmip6_ng_adm1_stagg")

# read in overview of model runs used
df_modelscen <- read_csv(file.path("data", "input", "df_modelscen.csv")) %>%
  mutate(identifier = paste0(model, "_", scenario, "_", ensemble))

# read in the GWL data
df_gwl <- read_csv(file.path("data", "input", "df_gwl.csv"))


################################################################################
########## PREPARE GLOBAL WARMING LEVELS DATA ##################################
################################################################################

# collapse start and finish years into a variable called window and clean up column names
df_gwl_table <- df_gwl %>% group_by(model, scenario, ensemble, warming_level) %>%
  summarise(beg = min(year), end = max(year)) %>%
  mutate(window = paste0(beg, "-", end)) %>%
  select(-beg, -end) %>%
  # convert to WIDE format (one column per GWL)
  pivot_wider(names_from = "warming_level", values_from = "window") %>%
  # value is NA for GWLs that a given model run does not reach - replace with an empty string
  mutate_at(vars("1", "1.5", "2", "3", "4"), ~ replace_na(.x, "")) %>%
  # add degC and + symbols to the warming levels in the column names
  setNames(names(.) %>% str_replace("(?<=[:digit:])$", "\u00B0C") %>%
             str_replace("^(?=[:digit:])", "\\+")) %>%
  # rename and rearrange columns
  select(GCM = "model", Scenario = "scenario", Ensemble = "ensemble", everything()) %>% arrange(GCM, Scenario)

# create a table with only one row per large ensemble
df_gwl_table %>%
  # discard the single run for CESM2 SSP3-7.0 (which is replaced with the entire large ensemble)
  filter(!(GCM == "CESM2" & Scenario == "ssp370")) %>%
  # ensure that we have only one row for the large ensembles
  filter(!(GCM == "CESM2-LE" & Ensemble != "r1i1p1f1")) %>%
  filter(!(GCM == "MPI-ESM1-2-LR" & Ensemble != "r1i1p1f1")) %>%
  # modify the information
  mutate(Ensemble = if_else(GCM %in% c("CESM2-LE", "MPI-ESM1-2-LR") & Scenario == "ssp370",
                            "Large ensemble", Ensemble)) %>%
  mutate_at(vars(contains("+")), ~ if_else(GCM %in% c("CESM2-LE", "MPI-ESM1-2-LR") & Scenario == "ssp370",
                                            "Varying", .x)) %>%
  # export to a .tex file
  xtable(label = "tab:gcm_ssp_list",
         digits = 0,
         caption = "Models, scenarios and ensemble members used for climatic projections") %>%
  print(file = file.path("..", "sharepoint", "Tables", paste0(Sys.Date(), " GCM_SSP_used_list.tex")),
        size="\\tiny")

# create a table only for the MPI large ensemble
df_gwl_table %>% 
  filter(GCM == "MPI-ESM1-2-LR" & Scenario == "ssp370") %>%
  # order by ensemble run (ascending)
  mutate(ensemble_run = str_extract(Ensemble, "(?<=^r)[:digit:]+(?=i)") %>% as.integer()) %>%
  arrange(ensemble_run) %>% select(-ensemble_run) %>%
  # export to a .tex file
  xtable(label = "tab:gcm_ssp_list",
         digits = 0,
         caption = "Models, scenarios and ensemble members used for climatic projections") %>%
  print(file = file.path("..", "sharepoint", "Tables", paste0(Sys.Date(), " GCM_SSP_used_list_mpi.tex")),
        size="\\tiny")

# create a table only for the MPI large ensemble
df_gwl_table %>% 
  filter(GCM == "CESM2-LE" & Scenario == "ssp370") %>%
  # order by ensemble run (ascending)
  mutate(ensemble_run = str_extract(Ensemble, "(?<=^r)[:digit:]+(?=i)") %>% as.integer()) %>%
  arrange(ensemble_run) %>% select(-ensemble_run) %>%
  # export to a .tex file
  xtable(label = "tab:gcm_ssp_list",
         digits = 0,
         caption = "Models, scenarios and ensemble members used for climatic projections") %>%
  print(file = file.path("..", "sharepoint", "Tables", paste0(Sys.Date(), " GCM_SSP_used_list_ensemble.tex")),
        size="\\tiny")


################################################################################
############### CALCULATE CEILING FOR BIAS-CORRECTED PRECIP DEVIATION ##########
################################################################################

# map all files for raw (= w/o bias correction) monthly precipitation totals
files_Wn_raw <- dir_climate_data_adm1 %>% str_replace("_stagg$", "_stagg_raw") %>% list.files(full.names = T, recursive = T) %>%
  str_subset("Wn") %>%
  # discard the single run for CESM2 SSP3-7.0 since we use the large ensemble
  str_subset("CESM2_ssp370", negate = T)

# check that we have 199 files as expected
if(length(files_Wn_raw) != 199) stop("Expected 199 files for raw monthly precipitation totals but got a different length. Please inspect")

# read in all files, add filename as a column and order by monthly precip deviation
Wn_raw <- map_dfr(files_Wn_raw, ~ read_feather(.x) %>% mutate(filename = .x)) %>%
  arrange(desc(Wn))

# calculate maximum monthly precip total for warming up to +3C
Wn_raw_max_by_model_3deg <- Wn_raw %>%
  
  # discard observations with a deviation of less than +2 since we're only interested in the maximum - just to reduce computation times
  filter(Wn > 2) %>% 
  
  # extract model, scenario and ensemble realization from the filename
  mutate(model = str_extract(filename, "(?<=stagg_raw\\/)[^_]+(?=_)"),
         scenario = str_extract(filename, "ssp[:digit:][:digit:][:digit:]"),
         ensemble = str_extract(filename, "r[:digit:]+i[:digit:]+p[:digit:]+f[:digit:]+")) %>%
  mutate(identifier = paste0(model, "_", scenario, "_", ensemble)) %>%
  
  # ensure that we only use model runs ultimately used in the paper (= featured in df_modelscen)
  filter(identifier %in% df_modelscen$identifier) %>%
  
  # discard the single run for CESM2 SSP3-7.0 (which is replaced with the entire large ensemble)
  filter(!(model == "CESM2" & scenario == "ssp370")) %>%
  
  # merge in GWLs (many-to-many since we have multiple GID_1 regions and the same model-scenario year can be in different GWL windows)
  left_join(df_gwl, by = c("model", "scenario", "ensemble", "year"), relationship = "many-to-many") %>%
  filter(warming_level <= 3) %>%
  
  # ensure that we treat low-emissions scenario (from CESM2) & SSP3-7.0 (from CESM2-LE) as the same model
  mutate(model = if_else(model == "CESM2-LE", "CESM2", model)) %>%
  
  # extract the highest value per model
  group_by(model) %>% slice_max(Wn, n = 1, with_ties = F) %>% ungroup() %>%
  
  # order models and add a rank variable
  arrange(desc(Wn)) %>% mutate(rank = 1:n()) %>% as_tibble() %>%
  
  # select variables of interest and add a degC symbol to warming_level
  select(rank, model, scenario, ensemble, year, warming_level, GID_1, Wn) %>%
  mutate(warming_level = paste0(warming_level, "\u00B0C")) %>%
  
  # rename columns for the final table
  setNames(names(.) %>% str_replace("^Wn$", "Maximum monthly precip. deviation")) %>%
  rename(GWL = "warming_level",
         Model = "model",
         Scenario = "scenario",
         Ensemble = "ensemble",
         Year = "year",
         Rank = "rank") %>%
  
  # convert year to a character, so it does not feature digits in the table
  mutate(Year = as.character(Year))

# inspect the result
Wn_raw_max_by_model_3deg

# since the MPI-ESM1-2-LR large ensemble produces stark outliers, we use the 2nd-highest value as ceiling for bias-corrected data
ceiling_Wn_bc <- Wn_raw_max_by_model_3deg %>% filter(Model != "MPI-ESM1-2-LR") %>%
                    pull(`Maximum monthly precip. deviation`) %>% max()

# inspect the result
ceiling_Wn_bc

# save out the ceiling for subsequent scripts
saveRDS(ceiling_Wn_bc, file.path("data", "ceiling_Wn_bc.rds"))

# write raw values out as a table
Wn_raw_max_by_model_3deg %>%
  select(-Rank) %>%
  xtable(label = "tab:monthlyprecip_maximum_bymodel",
         digits = 1,
         caption = "Maximum values for monthly precipitation deviation in raw CMIP6 outputs for global warming up to +3°C. The upper bound for the bias-corrected monthly precipitation deviation is based on the maximum value for EC-Earth3-Veg-LR (i.e., 11.6) because maximum values produced by MPI-ESM1-2-LR before bias correction represent clear outliers across all CMIP6 models under consideration.") %>%
  print(file = file.path("..", "sharepoint", "Tables", paste0(Sys.Date(), " Wn_maximum_by_model.tex")),
        size="\\small")

# read in bias-corrected values and check share of distribution affected
files_Wn_bc <- dir_climate_data_adm1 %>% list.files(full.names = T, recursive = T) %>%
  str_subset("Wn") %>%
  # discard the single run for CESM2 SSP3-7.0 since we use the large ensemble, and ERA5
  str_subset("CESM2_ssp370", negate = T) %>%
  str_subset("era5", negate = T)

# check that we have 199 files as expected
if(length(files_Wn_bc) != 199) stop("Expected 199 files for bias-corrected monthly precipitation totals but got a different length. Please inspect")

# read in all files, add filename as a column and order by monthly precip deviation
map_dfr(files_Wn_bc, ~ read_feather(.x)) %>%
  summarise(perc_cutoff = 1 - mean(Wn > ceiling_Wn_bc)) %>% pull(perc_cutoff)

# clean up
rm(Wn_raw, Wn_raw_max_by_model_3deg, ceiling_Wn_bc, files_Wn_raw, files_Wn_bc)
gc()


################################################################################
###################### PREPARE SSP DATA ########################################
################################################################################

# read in projected country-level SSP GDP
df_ssp_raw <- read_excel(file.path("data", "input", "220803 SSP-Database GDP PPP.xlsx"),
                     sheet = "data",
                     range = "A1:X921") %>%
  pivot_longer(cols = starts_with("2"), names_to = "year", values_to = "gdp") %>%
  # add a dummy indicating that data is projection, not historical
  mutate(is_historical = FALSE)

# read in the historic data starting in row 922
df_ssp_historic <- read_excel(file.path("data", "input", "220803 SSP-Database GDP PPP.xlsx"),
                       sheet = "data",
                       range = "A922:Q1103") %>%
              pivot_longer(cols = c(starts_with("2"), starts_with("1")), names_to = "year", values_to = "gdp") %>%
              rename(Scenario = "Scenario (History)") %>%
              # mark that data is historical
              mutate(is_historical = TRUE)

# combine the two
df_ssp <- bind_rows(df_ssp_raw, map_dfr(paste0("SSP", 1:5), ~ df_ssp_historic %>% mutate(Scenario = .x))) %>%
          # make year a numeric and harmonize column names
          rename(GID_0 = "Region", ssp = "Scenario") %>%
          mutate(year = as.integer(year)) %>%
          # use projections for 2010 (= slightly higher coverage) to avoid two values for the same year
          filter(!(year == 2010 & is_historical)) %>%
          # create a log version of GDP for log-linear interpolation
          mutate(gdp_ln = log(gdp))

# check that we have exactly one value per country, year and scenario
if((df_ssp %>% group_by(ssp, GID_0, year) %>% summarise(n = n()) %>% filter(n != 1) %>% nrow()) > 0) stop("There are year-country-SSP pairings with zero or > 1 values") 

# create the data set to be filled
# step 1: create a data frame with all possible combinations of SSPs, ADM0-level countries and years
df_ssp_out <- crossing(ssp = unique(df_ssp$ssp), GID_0 = unique(df_ssp$GID_0),
                       year = min(df_ssp$year):max(df_ssp$year)) %>%
  # merge in the SSP data
  left_join(df_ssp %>% select(ssp, GID_0, year, gdp, is_historical, gdp_ln),
            by = c("ssp", "GID_0", "year")) %>%
  # for each country and SSP, order by year and perform linear interpolation on the log-GDP via na.approx()
  group_by(ssp, GID_0) %>% arrange(ssp, GID_0, year) %>%
  mutate(gdp_ln_interpolated = zoo::na.approx(gdp_ln, maxgap = 9)) %>%
  # convert back to normal scale, then discard rows with missing GDP information after interpolation
  mutate(gdp_annual = exp(gdp_ln_interpolated)) %>%
  filter(!is.na(gdp_annual)) %>% ungroup() %>%
  # calculate weights
  group_by(ssp, year) %>%
  mutate(adm0_weight = gdp_annual/sum(gdp_annual)) %>%
  ungroup()

# check summary stats
summary(df_ssp_out)

# save this out as an RDS file
saveRDS(df_ssp_out %>% select(ssp, GID_0, year, gdp_annual, adm0_weight), file.path("data", "df_ssp.rds"))


################################################################################
##################### CREATE THE WORLD SHAPEFILE ###############################
################################################################################

# create a small-scale polygon of countries from rnaturalearthdata
world <- ne_countries(scale = "medium", returnclass = "sf") %>% select(sovereignt, type, admin, geounit, GID_0 = "iso_a3", geometry)

# read in our shapefile to compare
sf_gadm1_augmented <- st_read(file.path("data", "adm1_augmented.shp"))

# check for NA GID_0 codes
world %>% filter(is.na(GID_0))
# -> Northern Cyprus and Kosovo are classified as sovereign countries but have missing GID_0

# look up GID_0 codes in GADM
sf_gadm1_augmented %>% filter(str_detect(NAME_0, "Kosovo|Northern Cyprus"))
# -> XKO and XNC. We exclude the latter as the country is recognized only by Turkey

# add the codes
world$GID_0[world$sovereignt == "Kosovo"] <- "XKO"
#world$GID_0[world$sovereignt == "Northern Cyprus"] <- "XNC"

# inspect again
world %>% filter(str_detect(sovereignt, "Kosovo|Cyprus"))

# create a vector of all sovereign GID_0 territories
sovereign_countries_gid0 <- world %>% filter(type %in% c("Country", "Sovereign country") & admin == sovereignt & !is.na(GID_0)) %>% pull(GID_0)

# inspect tiny countries omitted in the shapefile pulled above
world_tinycountries <- ne_countries(type = "tiny_countries", returnclass = "sf", scale = "medium") %>%
  select(sovereignt, type, admin, GID_0 = "iso_a3", geounit, geometry)
world_tinycountries %>% filter(!GID_0 %in% world$GID_0)
# -> only Tuvalu (TUV) is categorized as a sovereign country

# add Tuvalu to the list of sovereign countries and the world shapefile
sovereign_countries_gid0 <- c(sovereign_countries_gid0, "TUV") %>% unique()
world <- bind_rows(world, world_tinycountries %>% filter(GID_0 == "TUV"))

# discard other small regions not included in our shapefile 
sovereign_countries_gid0 <- sovereign_countries_gid0[sovereign_countries_gid0 %in% sf_gadm1_augmented$GID_0]

# clean up the environment
rm(world_tinycountries)

# save out
saveRDS(world, file.path("data", "world.rds"))
saveRDS(sovereign_countries_gid0, file.path("data", "sovereign_countries_gid0.rds"))

