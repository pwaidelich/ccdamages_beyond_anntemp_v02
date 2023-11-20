rm(list = ls())

# load stagg (NOTE: package needs to be installed from Github since as of November 2023, it's not on CRAN yet: https://github.com/tcarleton/stagg)
library(stagg)
library(feather)

# load the ERA5 grid raster (included in the stagg package)
data(era5_grid)

# load the polygon shapefile for ADM1 regions
sf_gadm1_augmented <- st_read(file.path("data", "adm1_augmented.shp"))

# load the ERA5 grid raster (included in stagg)
data(era5_grid)

# calculate ADM1-level area weights
adm1_weights <- overlay_weights(
  
  polygons = sf_gadm1_augmented,
  
  polygon_id_col = "GID_1",
  
  grid = era5_grid,
  
  secondary_weights = NULL # no additional weighting aside from grid cell area
) 

# save out
write_feather(adm1_weights %>% as_tibble(), file.path("data", "df_adm1_weights_stagg.feather"))
