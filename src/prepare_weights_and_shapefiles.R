# clean the environment
rm(list = ls())

# load packages
library(tidyverse)
library(sf)         # for geodata operations


################################################################################
########## CREATE AUGMENTED ADM1-LEVEL DATA OBJECTS ############################
################################################################################

# The ADM1-level (geo-)data from GADM does not feature countries that have no ADM1-level subnational regions.
# Among countries in the SSP database, this is the case for the Maldives (MDV), Aruba (ABW) and Kiribati (KIR).
# Therefore, we create an augmented version of the ADM1-level data that features these countries with a single
# (imputed) ADM1-level region, which is assigned a GDP weight of one for ADM1-to-ADM0 aggregation.

# load ADM1 shapefile
sf_gadm1 <- st_read('../sharepoint/Data/GADM/gadm36_levels_gpkg/gadm36_levels.gpkg',
                    layer = "level1")

# load ADM0 shapefile
sf_gadm0 <- st_read('../sharepoint/Data/GADM/gadm36_levels_gpkg/gadm36_levels.gpkg',
                    layer = "level0")

# confirm that the three countries are not in the ADM1 shaefile
if(sum(str_count(sf_gadm1$GID_0, "MDV|ABW|KIR")) != 0) stop("There are pattern matches for MDV, ABW or KIR in the ADM1-level shapefile. Please inspect")

# create an augmented ADM0-level shapefile by adding these regions
sf_gadm1_augmented <- bind_rows(sf_gadm0 %>% filter(GID_0 %in% c("MDV", "ABW", "KIR")) %>% mutate(GID_1 = paste0(GID_0, ".1_1"),
                                                                                                  NAME_1 = paste0(NAME_0, " (imputed)")),
                                sf_gadm1)

# export the shapefile
st_write(sf_gadm1_augmented, file.path("data", "adm1_augmented.shp"), delete_layer = TRUE)



################################################################################
############################# ADD GDP WEIGHTS FOR ADDED REGIONS ################
################################################################################

# load the GDP weights created previously
gdp_weights <- read_csv(file.path("data", "adm1-weights.csv")) 

# add the three missing countries w/o ADM1-level subnational regions w/ a weight of 1
gdp_weights <- bind_rows(gdp_weights,
                         # add the countries with an imputed GDP of 1 (value is irrelevant as the region will have a weight of 100% regardless)
                         sf_gadm1_augmented %>% filter(GID_0 %in% c("ABW", "KIR", "MDV")) %>% as_tibble() %>% select(GID_0, GID_1, NAME_0, NAME_1) %>% mutate(GDP = 1)) %>% 
               # remove redundant rows (to avoid that by re-running the script, we add the countries multiple times)  
               distinct()

# throw an error if we do not get exactly 3 rows for these countries or if we do not have 3,613 rows
if(nrow(gdp_weights %>% filter(GID_0 %in% c("ABW", "KIR", "MDV"))) != 3) stop("Either more or less than 3 entries for ABW, KIR and MDV. Please inspect")
if(nrow(gdp_weights) != 3613) stop("gdp_weights should have 3,613 after adding the 3 extra regions. Please inspect")

# write out the augmented GDP weights file
write_csv(gdp_weights, file.path("data", "adm1-weights.csv"))



################################################################################
################## IDENTIFY FIRST YEAR OF GWL +0.84 ACROSS MODELS ##############
################################################################################

# NOTE: all projections are based only on years in the +0.84C GWL (= baseline) or
# subsequent years. Therefore, we identify the first year falling into this global
# warming level period across all models considered to facilitate cropping all
# climate data files automatically to the period starting from this year

# load GWL data
df_gwl <- read_csv(file.path("data", "input", "df_gwl.csv"))

# extract the very first global warming level year for the baseline period
first_year_GWL084 <- df_gwl %>% filter(warming_level == 0.84) %>% pull(year) %>% min()

# save out
saveRDS(first_year_GWL084, file.path("data", "first_year_GWL084.rds"))