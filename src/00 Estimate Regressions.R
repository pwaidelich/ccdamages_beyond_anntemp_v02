# NOTE: this script is forked from the 'main_results.R' in the Zenodo repo for the replication of Kotz et al, 2022
# several lines of code have been taken directly from this script to ensure that we fully replicate their regression
# results that inform our damage projections

# clean the environment
rm(list = ls())

# load required packages
library(tidyverse)
library(plm)
library(clubSandwich)
library(stargazer)
library(gsubfn)

# map the location of our climate indicators (for Tx5d)
dir_climate_data_adm1 <- file.path("/net", "cfc", "landclim1", "pwaidelich", "cmip6_ng_adm1_stagg")

# source helper functions
source("src/utils/analyse_impacts.R")
source("src/utils/project_impacts.R")
source("src/utils/extract_climate_indices.R")

#function to help stargazer package format in scientific format nicely, adapted from: https://stackoverflow.com/a/56924401
replace_numbers = function(x, low=0.01, high=1e3, digits = 3, scipen=-7, ...) {
  x = gsub(mark,'.',x)
  x.num = as.numeric(x)
  if (x.num<0.01){
    ifelse(
      (x.num >= low) & (x.num < high),
      round(x.num, digits = 2),
      prettyNum(x.num, digits=2, scientific = scipen, ...)
    )
  }
  else {
    ifelse(
      (x.num >= low) & (x.num < high),
      round(x.num, digits = digits),
      prettyNum(x.num, digits=digits, scientific = scipen, ...)
    )
  }
}



################################################################################
#################### REPRODUCE KOTZ ET AL 2022 MAIN SPEC #######################
################################################################################

#load data from Kotz et al 2022
table <- read.csv(file.path("data", "input", "PT_master.csv"))

# generate quadratic terms for annual total, monthly deviations, and the number of wet days
table$Pt2 <- table$Pt^2
table$Wn_2 <- table$Wn^2
table$wet_days_1_2 <- table$wet_days_1^2

# main econometric specification from Kotz et al., 2022
mod <- plm(dlgdp_pc_usd ~ Tstd + TmeanD + TmeanD:Tmean + TmeanLD + TmeanLD:TmeanL + Pt + Pt2 + Wn + Wn_2 + wet_days_1 + wet_days_1_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean,data=table,index=c('ID',"year"),model='within',effect='twoways')
coefs <- coef(mod)
cov_id <- vcovCR(mod,cluster=table$ID,type='CR0')
cov_iso <- vcovCR(mod,cluster=table$iso,type='CR0')
se_id <- sqrt(diag(cov_id))
se_iso <- sqrt(diag(cov_iso))

# clean up variable names and save out point estimates for coefficients
names(coefs) <- str_replace(names(coefs), ":", ".")
saveRDS(coefs, file = file.path("data", "kotzetal2022_coefs_main.rds"))

# clean up the clubSandwich objects into regular matrices
clean_up_clubsandwich <- function(clubsand_obj = NULL) {
  # set all attributes to NULL
  for(attr_name in c("type", "cluster", "bread", "v_scale",
                     "est_mats", "adjustments", "target", "inverse_var",
                     "ignore_FE")) attr(clubsand_obj, attr_name) <- NULL
  
  # set class to matrix and clean up rownames
  class(clubsand_obj) <- "matrix"
  rownames(clubsand_obj) <- str_replace(rownames(clubsand_obj), ":", ".")
  colnames(clubsand_obj) <- rownames(clubsand_obj)
  
  return(clubsand_obj)
}

# clean up the covariance matrices
cov_id <- clean_up_clubsandwich(cov_id)
cov_iso <- clean_up_clubsandwich(cov_iso)

# ensure that the order is consistent between the objects
if(!identical(rownames(cov_id), names(coefs))) stop("Rownames in cov_id do not match the names in the coefs vector")
if(!identical(rownames(cov_iso), names(coefs))) stop("Rownames in cov_iso do not match the names in the coefs vector")

# export as RDS files
saveRDS(cov_id, file = file.path("data", "kotzetal2022_cov_clustregional_main.rds"))
saveRDS(cov_iso, file = file.path("data", "kotzetal2022_cov_clustnational_main.rds"))



################################################################################
#################### PRODUCE KOTZ ET AL 2022 PERMUTATIONS ######################
################################################################################

# load our ERA5-based Tx5d values (area-weighted)
df_tx5d_era5 <- read_feather(file.path(dir_climate_data_adm1, "era5_084", "tx5d.feather"))

# merge Tx5d into table
table <- table %>% mutate(GID_1 = paste0(iso, ".", id_1, "_1")) %>% left_join(df_tx5d_era5, by = c("GID_1", "year"))

# add tx5d & pdsi from Callahan & Mankin 2022 data (ERA5, population-weighted)
df_tx5d_callahan <- read_csv(file.path("data", "input", "CallahanMankin2022_extremes_growth_panel_monthedd_1979-2016.csv")) %>%
  rename(year = "time", GID_1 = "region")

# merge Tx5d and PDSI from their data into table
table <- table %>% left_join(df_tx5d_callahan %>% select(GID_1, year, tx5d_callahan = "tx5d", pdsi_callahan = "pdsi"), by = c("year", "GID_1"))

# create a list storing the models, their coefficients, and their variance-covariance matrices
formula_list <- list()
mod_list <- list()
coef_list <- list()
cov_id_list <- list()
cov_iso_list <- list()

# create all the formulas of interest by iteratively adding regressors and name them accordingly
tempprecip_element <- "dlgdp_pc_usd ~ TmeanD + TmeanD:Tmean + TmeanLD + TmeanLD:TmeanL + Pt + Pt2"

# list all formulas to be used for regressions
formula_list <- list("T_Pt_Pt2zero" = as.formula(tempprecip_element),
                     "T_Pt_Tstd" = as.formula(paste0(tempprecip_element, " + Tstd")),
                     "T_Pt_Tstd_wetdays" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2")),
                     "T_Pt_Tstd_wetdays_Wn" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2 + Wn + Wn_2")),
                     "T_Pt_Tstd_wetdays_Wn_99p9" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2 + Wn + Wn_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean")),
                     "Kotz + tx5d" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2 + Wn + Wn_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean + tx5d + tx5d:Tmean")),
                     "Kotz + tx5d (linear)" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2 + Wn + Wn_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean + tx5d")),
                     "Kotz + tx5d (longterm)" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2 + Wn + Wn_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean + tx5d + tx5d:Tmean_longterm")),
                     "Kotz + pdsi (linear)" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2 + Wn + Wn_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean + pdsi_callahan")),
                     "Kotz + pdsi" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2 + Wn + Wn_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean + pdsi_callahan + pdsi_callahan:Tmean")),
                     "Kotz + pdsi (longterm)" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2 + Wn + Wn_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean + pdsi_callahan + pdsi_callahan:Tmean_longterm")),
                     "Kotz + tx5d_callahan" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2 + Wn + Wn_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean + tx5d_callahan + tx5d_callahan:Tmean"))
)

# estimate models by looping through the formula list
for(specname in names(formula_list)) {
  # NOTE: we add a long-term average of temperature for each region to interact Tx5d not only with current but also long-term temperatures
  mod_list[[specname]] <- plm(formula_list[[specname]],data=table %>% group_by(GID_1) %>% mutate(Tmean_longterm = mean(Tmean, na.rm = T)) %>% ungroup(),
                              index=c('ID',"year"),model='within',effect='twoways')
  
  # store coefficients and vcov matrices (for national- and regional clustering) into the list objects and add empty elements (= climate indicators that do not feature in the model)
  coef_list[[specname]] <- coef(mod_list[[specname]]) %>% setNames(names(.) %>% str_replace(":", "."))
  cov_id_list[[specname]] <- vcovCR(mod_list[[specname]],cluster=table$ID,type='CR0') %>% clean_up_clubsandwich()
  cov_iso_list[[specname]] <- vcovCR(mod_list[[specname]],cluster=table$iso,type='CR0') %>% clean_up_clubsandwich()
}

# check the sample size for all regression models
map_dbl(mod_list, .f = ~ length(.x$residuals))
# -> we lose some observations when using measures from Callahan & Mankin 2022

# save out the coefficient and the country-clustered variance-covariance matrix 
saveRDS(coef_list, file = file.path("data", "coef_list.rds"))
saveRDS(cov_id_list, file = file.path("data", "cov_id_list.rds"))
saveRDS(cov_iso_list, file = file.path("data", "cov_iso_list.rds"))


### Table with incrementally removing regressors from Kotz et al 2022 main spec

# export the respective regression tables
mod_for_table <- list(m1 = mod_list$T_Pt_Pt2zero,
                      m2 = mod_list$T_Pt_Tstd,
                      m3 = mod_list$T_Pt_Tstd_wetdays,
                      m4 = mod_list$T_Pt_Tstd_wetdays_Wn,
                      m5 = mod_list$T_Pt_Tstd_wetdays_Wn_99p9)
se_for_table <- list(m1 = sqrt(diag(cov_iso_list$T_Pt_Pt2zero)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                     m2 = sqrt(diag(cov_iso_list$T_Pt_Tstd)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                     m3 = sqrt(diag(cov_iso_list$T_Pt_Tstd_wetdays)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                     m4 = sqrt(diag(cov_iso_list$T_Pt_Tstd_wetdays_Wn)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                     m5 = sqrt(diag(cov_iso_list$T_Pt_Tstd_wetdays_Wn_99p9)) %>% setNames(names(.) %>% str_replace("\\.", ":")))

# print out the table
mark  = '::::'
star = stargazer(mod_for_table,se = se_for_table,dep.var.labels=c("Subnational income growth"),column.labels=NULL,omit.stat="F",align=TRUE,digits=10,decimal.mark  = mark, digit.separator = '',
                 order = c("^TmeanD$", "^TmeanLD$", "TmeanD\\:", "TmeanLD\\:", "Tstd", "^Pt$", "Pt2", "^Wn$", "Wn_2", "wet_days_1$", "wet_days_1_2",
                   "vwet_days1_am_99p9$", "Tmean\\:vwet_days1_am_99p9", "^tx5d$", "Tmean\\:tx5d", "tx5d\\:Tmean_longterm", "pdsi"),
                 label = "tab:regress_main",
                 title = "Regression models used to estimate dose-response functions")
cat(gsubfn(paste0("([0-9.\\-]+", mark, "[0-9.\\-]+)"), ~replace_numbers(x), star), sep='\n')

# clean up
rm(mod_for_table, se_for_table, star)


### Table for adding Tx5d (linear, interacted with temperature, etc.)

# export the respective regression tables
mod_for_table_tx5d <- list(m1 = mod_list$T_Pt_Tstd_wetdays_Wn_99p9,
                           m2 = mod_list$`Kotz + tx5d (linear)`,
                           m3 = mod_list$`Kotz + tx5d`,
                           m4 = mod_list$`Kotz + tx5d (longterm)`)

se_for_table_tx5d <- list(m1 = sqrt(diag(cov_iso_list$T_Pt_Tstd_wetdays_Wn_99p9)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                        m2 = sqrt(diag(cov_iso_list$`Kotz + tx5d (linear)`)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                     m3 = sqrt(diag(cov_iso_list$`Kotz + tx5d`)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                     m4 = sqrt(diag(cov_iso_list$`Kotz + tx5d (longterm)`)) %>% setNames(names(.) %>% str_replace("\\.", ":")))

# print out the table
star_tx5d = stargazer(mod_for_table_tx5d,se = se_for_table_tx5d,dep.var.labels=c("Subnational income per capita growth"),column.labels=NULL,omit.stat="F",align=TRUE,digits=10,decimal.mark  = mark, digit.separator = '',
                      order = c("^TmeanD$", "^TmeanLD$", "TmeanD\\:", "TmeanLD\\:", "Tstd", "^Pt$", "Pt2", "^Wn$", "Wn_2", "wet_days_1$", "wet_days_1_2",
                                "vwet_days1_am_99p9$", "Tmean\\:vwet_days1_am_99p9", "^tx5d$", "Tmean\\:tx5d", "tx5d\\:Tmean_longterm", "pdsi"),
                      title = "Including heat measured by Tx5d following Callahan & Mankin (2022)",
                      label = "tab:regress_tx5d")
cat(gsubfn(paste0("([0-9.\\-]+", mark, "[0-9.\\-]+)"), ~replace_numbers(x), star_tx5d), sep='\n')

# clean up
rm(mod_for_table_tx5d, se_for_table_tx5d, star_tx5d)


### Table for adding PDSI (linear, interacted with temperature, etc.)

# export the respective regression tables
mod_for_table_pdsi <- list(m1 = mod_list$T_Pt_Tstd_wetdays_Wn_99p9,
                           m2 = mod_list$`Kotz + pdsi (linear)`,
                           m3 = mod_list$`Kotz + pdsi`,
                           m4 = mod_list$`Kotz + pdsi (longterm)`,
                           m5 = mod_list$`Kotz + tx5d_callahan`)

se_for_table_pdsi <- list(m1 = sqrt(diag(cov_iso_list$T_Pt_Tstd_wetdays_Wn_99p9)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                          m2 = sqrt(diag(cov_iso_list$`Kotz + pdsi (linear)`)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                          m3 = sqrt(diag(cov_iso_list$`Kotz + pdsi`)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                          m4 = sqrt(diag(cov_iso_list$`Kotz + pdsi (longterm)`)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                          m5 = sqrt(diag(cov_iso_list$`Kotz + tx5d_callahan`)) %>% setNames(names(.) %>% str_replace("\\.", ":")))

# print out the table
star_pdsi = stargazer(mod_for_table_pdsi, se = se_for_table_pdsi, dep.var.labels=c("Subnational income per capita growth"),column.labels=NULL,omit.stat="F",align=TRUE,digits=10,decimal.mark  = mark, digit.separator = '',
                      order = c("^TmeanD$", "^TmeanLD$", "TmeanD\\:", "TmeanLD\\:", "Tstd", "^Pt$", "Pt2", "^Wn$", "Wn_2", "wet_days_1$", "wet_days_1_2",
                                "vwet_days1_am_99p9$", "Tmean\\:vwet_days1_am_99p9", "pdsi", "tx5d"),
                      title = "Including drought measured by the Palmer Drought Severity Index (PDSI)",
                      label = "tab:regress_pdsi")
cat(gsubfn(paste0("([0-9.\\-]+", mark, "[0-9.\\-]+)"), ~replace_numbers(x), star_pdsi), sep='\n')

# clean up
rm(mod_for_table_pdsi, se_for_table_pdsi, star_tx5d)

# NOTE: all regression tables are (slightly) reformatted manually, therefore, the
# the tables produced above reproduce the content of our SI tables but not exactly
# their final appearances

# calculate correlation coefficients for Tx5d
cor(table$tx5d, table$tx5d_callahan, use = "pair")
cor(table$Tmean, table$tx5d, use = "pair")



################################################################################
##################### KOTZ ET AL 2022 MAIN SPEC FOR INFLATION ADJUSTMENT #######
################################################################################

## a) use real GDPpc instead of nominal following Callahan & Mankin 2022

# load the GDP_deflator variable for the US and implement following Callahan & Mankin 2022
df_gdpdeflator <- read_csv(file.path("data", "input", "API_NY.GDP.DEFL.ZS_DS2_en_excel_v2_4353530.csv"))
df_deflator_us <- df_gdpdeflator %>% filter(`Country Code` == "USA") %>% select(-`Country Name`, -starts_with("Indicator")) %>%
                        gather(key = "year", value = "gdp_deflator", -c("Country Code")) %>% rename(GID_0 = "Country Code")

# extract 2010 value and calculate the factor by which we need to multiply to inflate/deflate to 2010USD
us_gdpdeflator_2010 <- df_deflator_us %>% filter(year == 2010) %>% pull(gdp_deflator)
df_deflator_us <- df_deflator_us %>% mutate(conversion_to_usd2010 = us_gdpdeflator_2010 / gdp_deflator,
                                            year = as.integer(year))

# clean up
rm(df_gdpdeflator, us_gdpdeflator_2010)

# quickly check that we can reproduce the dependent variable of Kotz et al 2022 (dlgdp_pc_usd) by taking ADM1-specific
# first-difference of log-transformed GDPpc
# (serves to rule out that data is not structured such that we cannot use lag() functions)
table %>% group_by(ID) %>% mutate(id_checking = 1:n(),
                                  lag_checking = if_else(id_checking == 1, NA_real_, dplyr::lag(lgdp_pc_usd)),
                                  dlgdp_pc_usd_checking = lgdp_pc_usd - lag_checking,
                                  diff_checking = dlgdp_pc_usd_checking - dlgdp_pc_usd) %>%
  select(id_checking, lgdp_pc_usd, lag_checking, dlgdp_pc_usd, dlgdp_pc_usd_checking, diff_checking) %>%
  summary()
# -> difference is basically zero (although not precisely due to precision issues)

# merge into Kotz et al 2022 data
table <- table %>% left_join(df_deflator_us %>% select(year, conversion_to_usd2010), by = "year") %>%
          # convert nominal GDPpc in current USD to 2010USD using the GDP deflator ratio for the US
          # NOTE: the Kotz et al 2022 data features some zero values for a Turkish region in gdp_pc_usd (with no implications for their analysis)
          #  to avoid numerical issues arising from taking log(0), we set these values to NA instead
          mutate(gdp_pc_2010 = if_else(gdp_pc_usd != 0, gdp_pc_usd*conversion_to_usd2010, NA_real_),
                 # create dependent variable by calculating logs
                 lgdp_pc_2010 = log(gdp_pc_2010)) %>%
          # take first-differences
          group_by(ID) %>% mutate(dlgdp_pc_2010 = lgdp_pc_2010 - dplyr::lag(lgdp_pc_2010, n = 1)) %>% ungroup()


# compare the two measures
table %>% select(dlgdp_pc_2010, dlgdp_pc_usd) %>% summary()
# -> very similar range but median and mean are considerably lower once we adjust for inflation

# repeat the main specification regression using inflation-adjuste GDPpc growth (in log-scale) as dep. var.
mod_infladjusted <- plm(dlgdp_pc_2010 ~ Tstd + TmeanD + TmeanD:Tmean + TmeanLD + TmeanLD:TmeanL + Pt + Pt2 + Wn + Wn_2 + wet_days_1 + wet_days_1_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean,data=table,index=c('ID',"year"),model='within',effect='twoways')

# extract variance-covariance matrix & st. errors for country-level clustering
cov_iso_infladjusted <- vcovCR(mod_infladjusted,cluster=table$iso,type='CR0')
se_iso_infladjusted <- sqrt(diag(cov_iso_infladjusted))

# ensure that order of variables is identical across coef()-output and the 'se_...' objects
if(!identical(names(coef(mod)), names(se_iso))) stop("order of climate indicators in coef(mod) and se_iso is not compatible")
if(!identical(names(coef(mod_infladjusted)), names(se_iso_infladjusted))) stop("order of climate indicators in coef(mod_infladjusted) and se_iso_infladjusted is not compatible")

# write out the inflation-adjusted regression table
star_infladjusted <- stargazer(mod, mod_infladjusted, se = list(se_iso, se_iso_infladjusted),dep.var.labels=c("Subnational income per capita growth"),column.labels=NULL,omit.stat="F",align=TRUE,digits=10,decimal.mark  = mark, digit.separator = '')
cat(gsubfn(paste0("([0-9.\\-]+", mark, "[0-9.\\-]+)"), ~replace_numbers(x), star_infladjusted), sep='\n')