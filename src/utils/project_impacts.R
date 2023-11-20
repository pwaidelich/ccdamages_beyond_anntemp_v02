# write a function to load the coefficients for the dose-response function
get_coefficients <- function(draw_from_vcov = FALSE, cluster_national = TRUE, seed = 2022,
                             n = 1, model_used = "main") {
  
  # check function arguments
  if(!is.logical(cluster_national)) stop("Argument cluster_national must be TRUE or FALSE")
  
  if(n != 1 & !draw_from_vcov) stop("Argument n can only be unequal to 1 if draw_from_vcov is TRUE")
  
  if(!model_used %in% c("main", "T_Pt_Pt2zero", "Kotz + tx5d", "T_Pt_Tstd")) {
    stop("model_used must be either main, T_Pt_Pt2zero, Kotz + tx5d, or T_Pt_Tstd")
  }
  
  # load the regression coefficients - if specification is not main, we load the list of coefficient vectors and extract the specified element by name (= model_used)
  if(model_used == "main") {
    coefs <- readRDS(file.path("data", "kotzetal2022_coefs_main.rds"))
  } else {
    coefs <- readRDS(file.path("data", "coef_list.rds"))[[model_used]]
  }
  
  # ensure that interactions are shown with .
  names(coefs) <- gsub(":", ".", names(coefs))
  
  # return the coefficients if no vcov-drawing is required
  if(!draw_from_vcov) return(coefs %>% as.matrix())
  
  # otherwise load the correct, draw the required sample and return it
  if(draw_from_vcov) {
    if(cluster_national) {
      if(model_used == "main") {
        vcov <- readRDS(file.path("data", "kotzetal2022_cov_clustnational_main.rds"))
      } else {
        vcov <- readRDS(file.path("data", "cov_iso_list.rds"))[[model_used]]
      }
    } else {
      if(model_used == "main") {
        vcov <- readRDS(file.path("data", "kotzetal2022_cov_clustregional_main.rds"))
      } else {
        vcov <- readRDS(file.path("data", "cov_id_list.rds"))[[model_used]]
      }
    }
    
    # round the matrix to avoid asymmetry due to precision, then draw
    vcov <- round(vcov, 18)
    set.seed(seed)
    out <- MASS::mvrnorm(n = n, mu = coefs %>% setNames(NULL), Sigma = vcov) %>% t()
    
    # NOTE: due to imprecision, mvrnorm() can return extremely small non-zero values for variables with coefficient & covariance set to zero
    # to address this, we manually overwrite these values with a clean zero
    vars_excluded <- names(coefs)[coefs == 0]
    out[vars_excluded, ] <- 0
    
    return(out)
  }
}


# function to calculate impacts in a given year relative to the baseline period
project_impacts <- function(data = NULL,
                            coefs = NULL,
                            coefs_Tonly = NULL,
                            df_gwl_baseline = NULL,
                            baseline_start_manual = 1979,
                            baseline_end_manual = 2019,
                            return_helper_columns = FALSE,
                            use_gwl_baseline = TRUE,
                            gwl_baseline_level = 0.84,
                            include_tx5d = FALSE,
                            return_as_datatable = F) {
  
  ## check inputs
  # names of the elements in the coefs object
  if(mean(c("Tstd", "TmeanD", "TmeanLD", "Pt", "Pt2", "Wn", "Wn_2", "wet_days_1", "wet_days_1_2","vwet_days1_am_99p9",
           "TmeanD.Tmean",  "TmeanLD.TmeanL", "Tmean.vwet_days1_am_99p9") %in% names(coefs)) != 1) stop("Names of elements in coefs do not match the required variables in Kotz et al 2022 regression")
  if(include_tx5d & mean(c("tx5d", "Tmean.tx5d") %in% names(coefs)) != 1) stop("tx5d & Tmean.tx5d must be elements in coefs if include_tx5d is TRUE")
  # column names in data
  if(mean(c("model", "scenario", "ensemble", "GID_1", "year", "Tmean", "Tstd", "Pt", "wet_days_1", "vwet_days1_am_99p9", "Wn") %in% names(data)) != 1) stop("Argument data must have columns named model, scenario, GID_1, year, Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn. Please inspect")
  if(include_tx5d & !"tx5d" %in% names(data)) stop("If include_tx5d is set to TRUE, data must feature a column called tx5d")
  # column names in df_gwl_baseline
  if(mean(c("model", "scenario", "warming_level", "ensemble", "year") %in% names(df_gwl_baseline)) != 1) stop("Argument df_gwl_baseline must feature columns called model, scenario, ensemble, warming_level, year")
  # argument types
  if(!"numeric" %in% class(coefs)) stop("Argument coefs must be of type numeric")
  if(!is.logical(return_helper_columns)) stop("Argument return_helper_columns must of type logical")
  # does data contain data for exactly one model-scenario pair?
  if((data %>% select(model, scenario) %>% distinct() %>% nrow()) != 1) stop("Argument data features more than one model-scenario pair. project_impacts() only works for data for a single model-scenario pair")
  # is the GWL baseline value a single numeric
  if(!gwl_baseline_level %in% df_gwl_baseline$warming_level) stop("Argument gwl_baseline_level must be a numeric that features in the warming_level column of df_gwl_baseline")
  
  # depending on user input, we extract baseline period definition either from user input or from the baseline global warming level window in df_gwl_baseline  
  if(!use_gwl_baseline) {
    
    baseline_start_used <- baseline_start_manual
    baseline_end_used <- baseline_end_manual
    
  } else {

    baseline_gwl_window <- data.table::as.data.table(data) %>%
        dtplyr::lazy_dt() %>%
        select(model, scenario, ensemble) %>% distinct() %>%
        left_join(df_gwl_baseline %>% filter(warming_level == gwl_baseline_level) %>% data.table::as.data.table(), by = c("model", "scenario", "ensemble")) %>%
        summarise(start = min(year), end = max(year)) %>%
        as_tibble()
    
    baseline_start_used <- baseline_gwl_window$start
    baseline_end_used <- baseline_gwl_window$end
    
  }
  
  if(!baseline_start_used %in% data$year | !baseline_end_used %in% data$year) stop("Argument data does not feature the specified baseline start/end years in column year")
  
  # identify baseline years via a dummy column
  data <- data.table::as.data.table(data) %>% mutate(is_baseline = year %in% baseline_start_used:baseline_end_used) %>%
    dtplyr::lazy_dt()
  
  # calculate the impacts in the baseline period
  df_base <- data %>%
    
    # subset to years in the baseline period
    filter(is_baseline) %>%
    
    # calculate the mean annual temperature across the baseline period (used to calculate extreme precip damages to avoid implicit adaptation through warming temperatures)
    group_by(GID_1) %>% mutate(Tmean_base = mean(Tmean)) %>% ungroup() %>%
    mutate(imp_Pt = Pt*coefs["Pt"] + (Pt)^2 * coefs["Pt2"],
           imp_Wn = Wn*coefs["Wn"] + (Wn)^2 * coefs["Wn_2"],
           imp_Tstd = Tstd * coefs["Tstd"],
           imp_wet_days_1 = wet_days_1 * coefs["wet_days_1"] + (wet_days_1)^2 * coefs["wet_days_1_2"],
           imp_vwet_days1_am_99p9 = vwet_days1_am_99p9 * coefs["vwet_days1_am_99p9"] + vwet_days1_am_99p9*Tmean_base*coefs["Tmean.vwet_days1_am_99p9"]
    ) %>%
    # calculate Tx5d impacts if required
    {if(include_tx5d) mutate(., imp_tx5d = tx5d * coefs["tx5d"] + tx5d*Tmean * coefs["Tmean.tx5d"]) else .}
  
  # take average in log-scale for the baseline period
  df_base_mean <- df_base %>%
    group_by(GID_1) %>% summarise_at(vars(starts_with("imp")), mean) %>%
    # add the Tmean_base (which is already averaged, see above)
    left_join(df_base %>% select(GID_1, Tmean_base) %>% distinct(), #%>% data.table::as.data.table(),
              by = "GID_1") #%>% as_tibble()
  
  # calculate impacts for all years
  df_out <- data %>%
    
    # discard all years prior to the baseline period as we do not use them anyways to save runtime
    filter(year >= baseline_start_used) %>%
    
    # select variables of interest
    select(model, scenario, ensemble, GID_1, year, is_baseline,
                                   Tmean, Tstd, Pt, wet_days_1, Wn, vwet_days1_am_99p9,
                                   any_of(c("tx5d"))) %>%
    # merge in the baseline average temperature (required to project out extreme precipitation damages)
    left_join(df_base_mean %>% select(Tmean_base, GID_1), # %>% data.table::as.data.table(),
              by = "GID_1") %>%
    
    # calculate current impacts (in log-scale)
    mutate(imp_Pt = Pt*coefs["Pt"] + (Pt)^2 * coefs["Pt2"],
           imp_Wn = Wn*coefs["Wn"] + (Wn)^2 * coefs["Wn_2"],
           imp_Tstd = Tstd * coefs["Tstd"],
           imp_wet_days_1 = wet_days_1 * coefs["wet_days_1"] + (wet_days_1)^2 * coefs["wet_days_1_2"],
           imp_vwet_days1_am_99p9 = vwet_days1_am_99p9 * coefs["vwet_days1_am_99p9"]  + vwet_days1_am_99p9 * Tmean_base * coefs["Tmean.vwet_days1_am_99p9"]
    ) %>%
    {if(include_tx5d) mutate(., imp_tx5d = tx5d * coefs["tx5d"] + tx5d*Tmean * coefs["Tmean.tx5d"]) else .} %>%
    
    # now merge in the average impacts during the baseline period (in log-scale) and denote them via a "_base" suffix
    left_join(df_base_mean %>% select(GID_1, starts_with("imp")), by = "GID_1", suffix = c("", "_base")) %>%
    
    # calculate differences (still in log-scale)
    mutate(imp_Tstd_diff = imp_Tstd - imp_Tstd_base,
           imp_Pt_diff = imp_Pt - imp_Pt_base,
           imp_Wn_diff = imp_Wn - imp_Wn_base,
           imp_wet_days_1_diff = imp_wet_days_1 - imp_wet_days_1_base,
           imp_vwet_days1_am_99p9_diff = imp_vwet_days1_am_99p9 - imp_vwet_days1_am_99p9_base
    ) %>%
    {if(include_tx5d) mutate(., imp_tx5d_diff = imp_tx5d - imp_tx5d_base) else .} %>%
    
    # temperature damages
    group_by(GID_1) %>% arrange(year, .by_group=T) %>%
    mutate(Tmean_fd = c(NA, diff(Tmean))) %>%
    mutate(dimp_temp = coefs["TmeanD"]*Tmean_fd + coefs["TmeanLD"]*lag(Tmean_fd) +
             coefs["TmeanD.Tmean"]*Tmean_fd*Tmean + coefs["TmeanLD.TmeanL"]*lag(Tmean_fd)*lag(Tmean)) %>%
    
    # impact is the cumulative sum for all years since the baseline start year
    mutate(imp_temp = cumsum(ifelse(is.na(dimp_temp), 0, dimp_temp))) %>%
    
    # impact over baseline is imp_temp minus the average over the baseline period
    mutate(imp_temp_during_base = if_else(is_baseline, imp_temp, NA_real_),
           imp_temp_diff = imp_temp - mean(imp_temp_during_base, na.rm = T)) %>%
    
    # repeat this for coefs_Tonly if they were provided by the user
    {if(!is.null(coefs_Tonly)) mutate(.,
                                      dimp_temp_Tonly = coefs_Tonly["TmeanD"]*Tmean_fd + coefs_Tonly["TmeanLD"]*lag(Tmean_fd) +
                                        coefs_Tonly["TmeanD.Tmean"]*Tmean_fd*Tmean + coefs_Tonly["TmeanLD.TmeanL"]*lag(Tmean_fd)*lag(Tmean)) %>%
        
        mutate(imp_temp_Tonly = cumsum(ifelse(is.na(dimp_temp_Tonly), 0, dimp_temp_Tonly))) %>%
        
        mutate(imp_temp_Tonly_during_base = if_else(is_baseline, imp_temp_Tonly, NA_real_),
               imp_temp_Tonly_diff = imp_temp_Tonly - mean(imp_temp_Tonly_during_base, na.rm = T)) else . } %>%
    
    # create total impacts and varextremes impact
    mutate(imp_total = imp_temp_diff + imp_Pt_diff + imp_Tstd_diff + imp_Wn_diff + imp_wet_days_1_diff + imp_vwet_days1_am_99p9_diff,
           imp_varextremes = imp_Tstd_diff + imp_Wn_diff + imp_wet_days_1_diff + imp_vwet_days1_am_99p9_diff) %>%
    
    # add Tx5d impacts to the totals if required
    {if(include_tx5d) mutate(., imp_total = imp_total + imp_tx5d_diff,
                                imp_varextremes = imp_varextremes + imp_tx5d_diff) else .} %>%
    
    # convert all impacts also from log-scale to %GDP
    mutate_at(vars(starts_with("imp")), ~ exp(.x) - 1) %>%
    
    # ungroup and convert back to a tibble depending on user input
    ungroup() %>%
    {if(!return_as_datatable) as_tibble(.) else .}
  
  # write out results (incl. helper columns if this is specified by the user)
  if(return_helper_columns) {
    
    return(df_out)
    
  } else {
    
    return(df_out %>% select(model, scenario, ensemble, GID_1, year,
                             imp_total, imp_varextremes, ends_with("_diff")))
  }
}


# wrapper to calculate impacts across Monte Carlo draws for the damage function
produce_gdp_projections <- function(model = NULL,
                                    scenario = NULL,
                                    ensemble = NULL,
                                    dir_climate_adm1 = NULL,
                                    dir_out_mc =  NULL, 
                                    overwrite_previous_files = FALSE,
                                    coefs_matrix_used = NULL,
                                    coefs_matrix_used_Tonly = NULL,
                                    adm1_weights_used = NULL,
                                    df_gwl_baseline = NULL,
                                    baseline_start_manual = 1979,
                                    baseline_end_manual = 2019,
                                    use_gwl_baseline = TRUE,
                                    gwl_baseline_level = 0.84,
                                    ceiling_precipdeviation = NULL,
                                    include_tx5d = FALSE,
                                    use_inner_parallelization = FALSE) {
  # input checks:
  if(!"matrix" %in% class(coefs_matrix_used)) stop("coefs_matrix_used must be of type matrix")
  if(!is.null(coefs_matrix_used_Tonly)) {
    if(ncol(coefs_matrix_used) != ncol(coefs_matrix_used_Tonly)) stop("coefs_matrix_used and coefs_matrix_used_Tonly must have the same no. of columns")
  }
  
  # create the identifier string
  identifier <- paste0(model, "_", scenario, "_", ensemble)
  
  # map the directory where ADM0 chunks should be saved and, if it's not there yet, create it
  dir_target_adm0 <- file.path(dir_out_mc, identifier)
  if(!dir.exists(dir_target_adm0)) dir.create(dir_target_adm0)
  
  # read in climate data
  df_adm1 <- collect_climdata_feathers(model = model,
                                       scenario = scenario,
                                       ensemble = ensemble,
                                       dir_files = dir_climate_adm1,
                                       include_tx5d = include_tx5d,
                                       simulate_Wn = F)
  
  # apply the ceiling to monthly precip deviation ('Wn') if specified 
  if(!is.null(ceiling_precipdeviation)) df_adm1 <- df_adm1 %>% mutate(Wn = if_else(Wn > ceiling_precipdeviation, ceiling_precipdeviation, Wn))
  
  # loop through the columns (= MC draws) in coefs_matrix_used via map()
  if(!use_inner_parallelization) {
    map(1:ncol(coefs_matrix_used), function(col_no) {
      
      # exit if the file already exists & overwriting is not required by the user
      # NOTE: we check for only one ADM0 code here since they're all written out simultaneously
      if(file.exists(file.path(dir_target_adm0,  paste0("ABW_mc", col_no, "_n", ncol(coefs_matrix_used), ".feather"))) & !overwrite_previous_files) {
        
        return(NULL)
      
      } else {
        # calculate impact projections using the col_no-th column of coefs_matrix_used
        project_impacts(data = df_adm1,
                        coefs = coefs_matrix_used[, col_no],
                        # NOTE: NULL[, i] returns NULL, so we do not need an ifelse() if coefs_matrix_used_Tonly is unspecified
                        coefs_Tonly = coefs_matrix_used_Tonly[, col_no],
                        df_gwl_baseline = df_gwl_baseline,
                        baseline_start_manual = baseline_start_manual,
                        baseline_end_manual = baseline_end_manual,
                        use_gwl_baseline = use_gwl_baseline,
                        gwl_baseline_level = gwl_baseline_level,
                        include_tx5d = include_tx5d,
                        return_helper_columns = F,
                        return_as_datatable = T) %>%
          
          # aggregate to ADM0 - this implicitly converts from data.table to tibble
          aggregate_adm1_to_adm0(adm1_weights_used = adm1_weights_used, is_mc_output = F, is_input_datatable = T) %>%
          
          # add the monte carlo draw no. as a column and reorder columns
          mutate(monte_carlo_draw = col_no) %>%
          select(model, scenario, ensemble, monte_carlo_draw, GID_0, year, everything()) %>%
          
          # prepare chunked out export by inserting the filename into the tibble
          mutate(filename = file.path(dir_target_adm0, paste0(GID_0, "_mc", col_no, "_n", ncol(coefs_matrix_used), ".feather"))) %>%
          
          # write out in chunks using group_walk (= basically a map by each group of the tibble)
          group_by(GID_0) %>%
          group_walk(~ .x %>% dplyr::select(-filename) %>% write_feather(path = .x$filename[1]), .keep = TRUE)
        
        return(TRUE)
      }
    })
  } else {
    # if inner parallelization is activated, we use a future_map instead of a map() to boost performance - NOTE: this is deprecated if we parallelize at a higher level (e.g. model-scenario)
    future_map(1:ncol(coefs_matrix_used), function(col_no) {
      
      # exit if the file already exists & overwriting is not required by the user
      if(file.exists(file.path(dir_target_adm0,  paste0("ABW_mc", col_no, "_n", ncol(coefs_matrix_used), ".feather"))) & !overwrite_previous_files) {
        
        return(NULL)
      
      } else {
        # calculate impact projections using the col_no-th column of coefs_matrix_used
        project_impacts(data = df_adm1,
                        coefs = coefs_matrix_used[, col_no],
                        # NOTE: NULL[, i] returns NULL, so we do not need an ifelse() if coefs_matrix_used_Tonly is unspecified
                        coefs_Tonly = coefs_matrix_used_Tonly[, col_no], 
                        df_gwl_baseline = df_gwl_baseline,
                        baseline_start_manual = baseline_start_manual,
                        baseline_end_manual = baseline_end_manual,
                        use_gwl_baseline = use_gwl_baseline,
                        gwl_baseline_level = gwl_baseline_level,
                        include_tx5d = include_tx5d,
                        return_helper_columns = F,
                        return_as_datatable = T) %>%
          
          # aggregate to ADM0 - this implicitly converts from data.table to tibble
          aggregate_adm1_to_adm0(adm1_weights_used = adm1_weights_used, is_mc_output = F, is_input_datatable = T) %>%
          
          # add the monte carlo draw no. as a column and reorder columns
          mutate(monte_carlo_draw = col_no) %>%
          select(model, scenario, ensemble, monte_carlo_draw, GID_0, year, everything()) %>%
          
          # prepare chunked out export by inserting the filename into the tibble
          mutate(filename = file.path(dir_target_adm0, paste0(GID_0, "_mc", col_no, "_n", ncol(coefs_matrix_used), ".feather"))) %>%
          
          # write out in chunks using group_walk (= basically a map by each group of the tibble)
          group_by(GID_0) %>%
          group_walk(~ .x %>% dplyr::select(-filename) %>% write_feather(path = .x$filename[1]), .keep = TRUE)
        
        return(TRUE)
      }
    })
  }
  
  return(TRUE)
}


# function to aggregate ADM1-level results to ADM0 using user-defined weights
aggregate_adm1_to_adm0 <- function(data = NULL, adm1_weights_used = NULL, is_mc_output = NULL, is_input_datatable = F) {
  
  if(!is_input_datatable) {
    # if input is a data.table object, meaning that reducing runtime is prioritized, we skip some input checks and hence only carry them out here
    if(mean(unique(data$GID_1) %in% adm1_weights_used$GID_1) < 1) stop("Not all ADM1 regions in data are covered by the ADM1 weights data provided")
    if(mean(c("GID_0", "GID_1", "adm1_weight") %in% names(adm1_weights_used)) != 1) stop("adm1_weights_used must contain columns GID_0, GID_1 and adm1_weight")
    
    # convert to data.table and initiate lazy_dt()
    data <- data.table::as.data.table(data) %>% dtplyr::lazy_dt()
  }
  
  if(is.null(is_mc_output)) stop("is_mc_output must be set to TRUE or FALSE depending on whether argument data is created via Monte Carlo analysis or not")
  
  df_out <- data %>%
            left_join(adm1_weights_used %>% select(GID_1, GID_0, adm1_weight) %>% data.table::as.data.table(), by = "GID_1") %>%
           {if(is_mc_output) group_by(., model, scenario, ensemble, monte_carlo_draw, year, GID_0) else group_by(., model, scenario, ensemble, year, GID_0)} %>%
           {if(!is_input_datatable) mutate(., test_dummy = 1) %>% summarise_at(vars(starts_with("imp"), "test_dummy"), ~ sum(.x * adm1_weight)) else summarise_at(., vars(starts_with("imp")), ~ sum(.x * adm1_weight))} %>%
            ungroup() %>%
            # collect final results
            as_tibble()
  
  # the aggregated test_dummy reflects to what weights sum up - throw an error if this is not one (rounding here to avoid precision issues)
  # again, we do this only if input is not a datatable
  if(!is_input_datatable) {
    if(mean(round(df_out$test_dummy, 9) == 1) != 1) stop(paste0("For some ADM0-level regions, weights sum up to the following values:\n",
                                                               paste0(df_out %>% filter(round(test_dummy, 10) != 1) %>% pull(test_dummy) %>% unique(), collapse = ", "),
                                                               "\nThis applies to model ", unique(df_out$model), " scenario ", unique(df_out$scenario)))
    
    return(df_out %>% select(-test_dummy))
  } else {
    return(df_out)
  }
}


# function to aggregate ADM1-level results to ADM0 using user-defined weights
aggregate_adm0_to_global <- function(data = NULL, adm0_weights_used = NULL, is_mc_output = NULL) {
  
  if(mean(c("GID_0", "ssp", "year", "adm0_weight") %in% names(adm0_weights_used)) != 1) stop("adm1_weights_used must contain columns GID_0, ssp, year, and adm0_weight")
  
  df_out <- data %>% mutate(ssp = paste0("SSP", str_extract(scenario, "(?<=^ssp)[:digit:]"))) %>%
            data.table::as.data.table() %>%
            dtplyr::lazy_dt() %>%
            # merge in GDP weights for the respective SSP
            left_join(adm0_weights_used %>% select(ssp, GID_0, year, adm0_weight) %>%
                        data.table::as.data.table(), by = c("ssp", "GID_0", "year")) %>%
            # exclude ADM0-level territories for which there is no GDP data in SSP Database and hence no weight
            filter(!is.na(adm0_weight)) %>%
            # group by model-scenario-year as well as MC run and calculate weighted averages (NOTE: data is in %GDP, not in log-scale)
            {if(is_mc_output) group_by(., model, scenario, ensemble, year, monte_carlo_draw) else group_by(., model, scenario, ensemble, year)} %>%
            mutate(test_dummy = 1) %>%
            summarise_at(vars(starts_with("imp"), "test_dummy"), ~ Hmisc::wtd.mean(.x, weights = adm0_weight, normwt = T)) %>% ungroup() %>%
            # collect final results
            as_tibble()
  
  # the aggregated test_dummy reflects to what weights sum up - throw an error if this is not one (rounding here to avoid precision issues)
  if(mean(round(df_out$test_dummy, 5) ==1) != 1) stop(paste0("For some years, weights sum up to the following values:\n",
                                                             paste0(df_out %>% filter(round(test_dummy, 10) != 1) %>% select(year, test_dummy) %>% distinct(), collapse = ", "),
                                                             "\nThis applies to model ", unique(df_out$model), " scenario ", unique(df_out$scenario)))
            
  # return the data frame without the test dummy variable
  return(df_out %>% select(-test_dummy))
}


# write a wrapper around aggregate_adm0_to_global() to facilitate parallelizing the ADM0-to-global aggregation
aggregate_to_global_parallel <- function(filename = NULL,
                                         dir_adm0 = NULL,
                                         dir_global = NULL,
                                         adm0_weights_used = NULL,
                                         overwrite_previous_files = FALSE) {
  
  # skip if the file already exists by returning NULL and not doing anything else
  if(file.exists(file.path(dir_global, filename)) & !overwrite_previous_files) {
    
    return(NULL)
    
    # else, read in the file, aggregate using the specified weights and saving out with the same filename in dir_global
  } else {
    
    read_feather(file.path(dir_adm0, filename)) %>%
      aggregate_adm0_to_global(adm0_weights_used = adm0_weights_used, is_mc_output = TRUE) %>%
      write_feather(file.path(dir_global, filename))
    
    return(TRUE)
  }
}



### OLD VERSION FROM THE INITIAL NCLIM SUBMISSION

# NOTE: we keep this simply to calculate the difference between averaging baseline period impacts
# in %GDP (= initial submission) or in log scale (= revision #1 version)

# function to calculate impacts in a given year relative to the baseline period
project_impacts_initial_submission <- function(data = NULL, coefs = NULL, df_gwl_baseline = NULL, baseline_start_manual = 1979, baseline_end_manual = 2019,
                            return_helper_columns = FALSE, use_gwl_baseline = TRUE, gwl_baseline_level = 0.84) {
  ## check inputs
  # names of the elements in the coefs object
  if(mean(c("Tstd", "TmeanD", "TmeanLD", "Pt", "Pt2", "Wn", "Wn_2", "wet_days_1", "wet_days_1_2","vwet_days1_am_99p9",
            "TmeanD.Tmean",  "TmeanLD.TmeanL", "Tmean.vwet_days1_am_99p9") %in% names(coefs)) != 1) stop("Names of elements in coefs do not match the required variables in Kotz et al 2022 regression")
  # column names in data
  if(mean(c("model", "scenario", "GID_1", "year", "Tmean", "Tstd", "Pt", "wet_days_1", "vwet_days1_am_99p9", "Wn") %in% names(data)) != 1) stop("Argument data must have columns named model, scenario, GID_1, year, Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn. Please inspect")
  # column names in df_gwl_baseline
  if(mean(c("model", "scenario", "warming_level", "year") %in% names(df_gwl_baseline)) != 1) stop("Argument df_gwl_baseline must feature columns called model, scenario, warming_level, year")
  # argument types
  if(!"numeric" %in% class(coefs)) stop("Argument coefs must be of type numeric")
  if(!is.logical(return_helper_columns)) stop("Argument return_helper_columns must of type logical")
  # does data contain data for exactly one model-scenario pair?
  if((data %>% select(model, scenario) %>% distinct() %>% nrow()) != 1) stop("Argument data features more than one model-scenario pair. project_impacts() only works for data for a single model-scenario pair")
  # is the GWL baseline value a single numeric
  if(!gwl_baseline_level %in% df_gwl_baseline$warming_level) stop("Argument gwl_baseline_level must be a numeric that features in the warming_level column of df_gwl_baseline")
  
  # depending on user input, we extract baseline period definition either from user input or from the baseline global warming level window in df_gwl_baseline  
  if(!use_gwl_baseline) {
    
    baseline_start_used <- baseline_start_manual
    baseline_end_used <- baseline_end_manual
    
  } else {
    
    baseline_gwl_window <- data %>% select(model, scenario) %>% distinct() %>%
      left_join(df_gwl_baseline %>% filter(warming_level == gwl_baseline_level), by = c("model", "scenario")) %>%
      summarise(start = min(year), end = max(year))
    
    baseline_start_used <- baseline_gwl_window$start
    baseline_end_used <- baseline_gwl_window$end
    
  }
  
  if(!baseline_start_used %in% data$year | !baseline_end_used %in% data$year) stop("Argument data does not feature the specified baseline start/end years in column year")
  
  # identify baseline years via a dummy column
  data <- data %>% mutate(is_baseline = year %in% baseline_start_used:baseline_end_used)
  
  # calculate the impacts in the baseline period
  df_base <- data %>% filter(is_baseline) %>%
    # calculate the mean annual temperature across the baseline period (used to calculate extreme precip damages to avoid implicit adaptation through warming temperatures)
    group_by(GID_1) %>% mutate(Tmean_base = mean(Tmean)) %>% ungroup() %>%
    mutate(imp_Pt = Pt*coefs["Pt"] + (Pt)^2 * coefs["Pt2"],
           imp_Wn = Wn*coefs["Wn"] + (Wn)^2 * coefs["Wn_2"],
           imp_Tstd = Tstd * coefs["Tstd"],
           imp_wet_days_1 = wet_days_1 * coefs["wet_days_1"] + (wet_days_1)^2 * coefs["wet_days_1_2"],
           imp_vwet_days1_am_99p9 = vwet_days1_am_99p9 * coefs["vwet_days1_am_99p9"] + vwet_days1_am_99p9*Tmean_base*coefs["Tmean.vwet_days1_am_99p9"]
    )
  
  # take average (NOTE: we convert to %GDP scale for this to avoid adding implicit risk premia by averaging the log impacts)
  df_base_mean <- df_base %>% group_by(GID_1) %>% summarise_at(vars(starts_with("imp")), ~ log(mean(exp(.x)))) %>%
    # add the Tmean_base (which is already averaged, see above)
    left_join(df_base %>% select(GID_1, Tmean_base) %>% distinct(), by = "GID_1")
  
  # calculate impacts for all years
  df_out <- data %>%
    # discard all years prior to the baseline period as we do not use them anyways to save runtime
    filter(year >= baseline_start_used) %>%
    # select variables of interest
    select(model, scenario, GID_1, year, is_baseline,
           ## FOLLOWING LINE WAS INSERTED POST INITIAL SUBMISSION TO ENSURE CONSISTENCY WITH CHANGES TO extract_climate_indices.R
           ensemble,
           ## END OF ADDED LINE
           Tmean, Tstd, Pt, wet_days_1, Wn, vwet_days1_am_99p9) %>%
    # merge in the baseline average temperature (required to project out extreme precipitation damages)
    left_join(df_base_mean %>% select(Tmean_base, GID_1), by = "GID_1") %>%
    # calculate current impacts (in log-scale)
    mutate(imp_Pt = Pt*coefs["Pt"] + (Pt)^2 * coefs["Pt2"],
           imp_Wn = Wn*coefs["Wn"] + (Wn)^2 * coefs["Wn_2"],
           imp_Tstd = Tstd * coefs["Tstd"],
           imp_wet_days_1 = wet_days_1 * coefs["wet_days_1"] + (wet_days_1)^2 * coefs["wet_days_1_2"],
           imp_vwet_days1_am_99p9 = vwet_days1_am_99p9 * coefs["vwet_days1_am_99p9"]  + vwet_days1_am_99p9 * Tmean_base * coefs["Tmean.vwet_days1_am_99p9"]
    ) %>%
    # now merge in the average impacts during the baseline period (in log-scale) and denote them via a "_base" suffix
    left_join(df_base_mean %>% select(GID_1, starts_with("imp")), by = "GID_1", suffix = c("", "_base")) %>%
    # calculate differences (still in log-scale)
    mutate(imp_Tstd_diff = imp_Tstd - imp_Tstd_base,
           imp_Pt_diff = imp_Pt - imp_Pt_base,
           imp_Wn_diff = imp_Wn - imp_Wn_base,
           imp_wet_days_1_diff = imp_wet_days_1 - imp_wet_days_1_base,
           imp_vwet_days1_am_99p9_diff = imp_vwet_days1_am_99p9 - imp_vwet_days1_am_99p9_base
    ) %>%
    # temperature damages
    group_by(GID_1) %>% arrange(year, .by_group=T) %>%
    mutate(Tmean_fd = c(NA, diff(Tmean))) %>%
    mutate(dimp_temp = coefs["TmeanD"]*Tmean_fd + coefs["TmeanLD"]*lag(Tmean_fd) +
             coefs["TmeanD.Tmean"]*Tmean_fd*Tmean + coefs["TmeanLD.TmeanL"]*lag(Tmean_fd)*lag(Tmean)) %>%
    # impact is the cumulative sum for all years since the baseline start year
    mutate(imp_temp = cumsum(ifelse(is.na(dimp_temp), 0, dimp_temp))) %>%
    # impact over baseline is imp_temp minus the average over the baseline period (in absolute scale, not log scale)
    mutate(imp_temp_absolutescale_base = if_else(is_baseline, exp(imp_temp), NA_real_),
           imp_temp_diff = imp_temp - log(mean(imp_temp_absolutescale_base, na.rm = T))) %>%
    # create total impacts and varextremes impact
    mutate(imp_total = imp_temp_diff + imp_Pt_diff + imp_Tstd_diff + imp_Wn_diff + imp_wet_days_1_diff + imp_vwet_days1_am_99p9_diff,
           imp_varextremes = imp_Tstd_diff + imp_Wn_diff + imp_wet_days_1_diff + imp_vwet_days1_am_99p9_diff) %>%
    # convert all impacts also from log-scale to %GDP, so we can average and aggregate them up
    mutate_at(vars(starts_with("imp")), ~ exp(.x) - 1) %>%
    ungroup()
  
  # write out results (incl. helper columns if this is specified by the user)
  if(return_helper_columns) {
    return(df_out)
  } else {
    return(df_out %>% select(model, scenario, GID_1, year,
                             imp_total, imp_varextremes, imp_temp_diff, imp_Tstd_diff, imp_Pt_diff, imp_Wn_diff,
                             imp_wet_days_1_diff, imp_vwet_days1_am_99p9_diff))
  }
}
