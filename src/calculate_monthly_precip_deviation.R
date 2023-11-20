rm(list = ls())

library(dtplyr)
library(feather)
library(tidyverse)

# map ADM1-level climate information
dir_climate_data_adm1 <- file.path("/net", "cfc", "landclim1", "pwaidelich", "cmip6_ng_adm1_stagg")

# load GWL data (since baseline of monthly precip deviations is GWL-dependent for CMIP6 models)
df_gwl <- read_csv(file.path("data", "input", "df_gwl.csv"))

# load overview of model runs considered
df_modelscen <- read_csv(file.path("data", "input", "df_modelscen.csv")) %>%
  # combine model, scenario and ensemble into an identifier used for folder names
  mutate(identifier = paste0(model, "_", scenario, "_", ensemble))

# write the function to calculate standardized monthly precip deviation following Kotz et al Methods section
calculate_monthly_precip_deviation <- function(filepath = NULL, # filepath of .feather w/ monthly precipitation totals 
                                               is_era5 = FALSE, # dummy indicating that .feather is ERA5, not CMIP6
                                               df_gwl_baseline = NULL, # tibble with GWL (= use df_gwl)
                                               gwl_baseline_level = 0.84 # GWL period used as baseline for monthly precip totals
                                               ) {
  
  # abort if filepath does not exist
  if(!file.exists(filepath)) stop("filepath does not exist")
  
  # skip if the file to be created already exists
  if(file.exists(filepath %>% str_replace("mon_totd\\.feather", "Wn.feather"))) {
    
    paste0(filepath, " already exists. Skipping...") %>% print()
    return(NULL)
    
  }
  
  # identify model, scenario and ensemble based on the filepath if we do not deal with ERA5
  if(!is_era5) {
    
    identifier <- str_extract(filepath, "[^\\/]+(?=\\/mon_totd\\.feather$)")
    model_used <- str_extract(identifier, "^[^_]+")
    scenario_used <- str_extract(identifier, "ssp[:digit:][:digit:][:digit:]")
    ensemble_used <- str_extract(identifier, "(?<=ssp[:digit:][:digit:][:digit:]_)r[:digit:]+i[:digit:]+p[:digit:]+f[:digit:]+")
  }
  
  # read in and separate the year variable (which has format "2019-01" for Jan 2019) into month and year
  df <- read_feather(filepath) %>%
    mutate(month = str_extract(year, "(?<=-)[01][:digit:]$"),
           year = str_extract(year, "^[12][:digit:][:digit:][:digit:]"))
  
  # throw an error if we find more/less than 12 months in the data
  if(length(unique(df$month)) != 12) stop(paste0("Less (or more) than 12 months detected in ", filepath, ". Aborting"))
  
  # NOTE: for ERA5 1979-2019 files, the entire file is used as baseline - otherwise we subset to the specified baseline GWL
  if(!is_era5) {
    
    # subset the GWL data to the baseline period for this model run
    df_gwl_baseperiod <- df_gwl_baseline %>% filter(model == model_used, scenario == scenario_used, ensemble == ensemble_used, warming_level == gwl_baseline_level)
    
    # throw an error if we do not have exactly 41 years - otherwise subset df accordingly, convert to data.table & initiate lazy_dt()
    if(nrow(df_gwl_baseperiod) != 41) stop("df_gwl_baseperiod must feature exactly 41 years")
    df_base <- data.table::as.data.table(df) %>% dtplyr::lazy_dt() %>% filter(year %in% df_gwl_baseperiod$year)
  
  }  else {
    
    # for ERA5, simply convert to data.table and initiate lazy_dt()
    df_base <- data.table::as.data.table(df) %>% dtplyr::lazy_dt()
  }
  
  # calculate mean & st.dv. by month & region across baseline period
  df_monthbase <- df_base %>%
    group_by(GID_1, month) %>%
    summarise(month_average_base = mean(mon_totd),
              month_stdv_base = sd(mon_totd)) %>% ungroup()
  
  # calculate annual precipitation by year & region, then average across baseline
  df_annualbase <- df_base %>%
    group_by(GID_1, year) %>%
    summarise(ann_totd = sum(mon_totd)) %>% ungroup() %>%
    group_by(GID_1) %>% summarise(annual_average_base = mean(ann_totd)) %>%
    ungroup()
  
  # calculate monthly precip deviation
  df_out <- df %>%
    
    data.table::as.data.table() %>%
    lazy_dt() %>%
    
    # merge in baseline values calculated above
    left_join(df_monthbase, by = c("GID_1", "month")) %>%
    left_join(df_annualbase, by = "GID_1") %>%
    
    # assign a zero value for months w/o precipitation (to avoid numerical issues with dividing by st.dv.)
    # since these months receive a weight of zero through month_average_base/annual_average_base
    mutate(Wn_month = if_else(month_stdv_base == 0, 0,
                              ((mon_totd - month_average_base)/month_stdv_base)*(month_average_base/annual_average_base))) %>%
    
    # add up for all months within the same year
    group_by(GID_1, year) %>% summarise(Wn = sum(Wn_month)) %>%
    
    # collect results of data.table
    as_tibble() %>%
    
    # convert year from a character to an integer
    mutate(year = as.integer(year))
  
  # write out into the same directory
  write_feather(df_out, filepath %>% str_replace("mon_totd\\.feather", "Wn.feather"))
  
  return(TRUE)
}

# apply to ERA5 monthly precipitation totals  
calculate_monthly_precip_deviation(filepath = file.path(dir_climate_data_adm1, "era5_084", "mon_totd.feather"), is_era5 = T)

# apply to all CMIP6 models for which ADM1-level monthly precipitation totals are available
map(1:nrow(df_modelscen),
    .f = function(row_number) {

      # print status to track progress
      print(row_number)
      paste0("Calculating Wn for ", df_modelscen$identifier[row_number], "...") %>% print()

      # map the filepath based on the identifier column in df_modelscen
      filepath = file.path(dir_climate_data_adm1, df_modelscen$identifier[row_number], "mon_totd.feather")

      # calculate monthly precipitation deviation if monthly precip totals .feather file exists
      if(file.exists(filepath)) {

       calculate_monthly_precip_deviation(filepath = filepath, is_era5 = F, df_gwl_baseline = df_gwl)

      } else {
        
        print("File already exists. Skipping...")
        return(NULL)
      }

  })


# repeat for raw model outputs (to calibrate the upper bound on bias-corrected monthly precip deviation)
map(1:nrow(df_modelscen),
    .f = function(row_number) {

      # print status to track progress
      print(row_number)
      paste0("Calculating Wn w/o bias-correction for ", df_modelscen$identifier[row_number], "...") %>% print()

      # map the filepath - NOTE: now we use climate data in the "..._stagg_raw" folder (= w/o bias correction)
      filepath = file.path(str_replace(dir_climate_data_adm1, "_stagg$", "_stagg_raw"), df_modelscen$identifier[row_number], "mon_totd.feather")

      if(file.exists(filepath)) {

        calculate_monthly_precip_deviation(filepath = filepath, is_era5 = F, df_gwl_baseline = df_gwl)

      } else {
        print("File already exists. Skipping...")
        return(NULL)
      }

})

# apply to the single MPI-ESM1-2-LR model run used for robustness check with alternative bias correction
calculate_monthly_precip_deviation(file.path(dir_climate_data_adm1 %>% str_replace("_stagg", "_stagg_altbiascor") ,
                                             "MPI-ESM1-2-LR_ssp370_r1i1p1f1", "mon_totd.feather"),
                                   df_gwl_baseline = df_gwl,
                                   is_era5 = F) 
