# Climate damage projections beyond annual temperature
Replication materials for Waidelich, Batibeniz, Rising, Kikstra, and Seneviratne (2023). If you have questions about the code or find any errors/bugs in it, please reach out to Paul Waidelich at paul.waidelich(at)gess.ethz.ch (corresponding author).

NOTE: THIS CLEAN REPOSITORY SERVES FOR THE REVIEW PROCESS ONLY TO ALLOW FOR TRANSPARENCY REGARDING THE SCRIPTS & ANALYSES UNDERLYING THE SUBMITTED MANUSCRIPT. AN EXTENDED REPOSITORY THAT ALSO COVERS LARGE INPUT & OUTPUT DATA FILES (STORED ON THE `sharepoint` FOLDER AND THE EXTERNAL DRIVE MENTIONED IN THE NEXT SECTION OF THIS README) WILL BE PUBLISHED UPON ACCEPTANCE.

## Organization of the overall project repository
Please note that both the underlying climate data and the subnational GDP projections produced require a substantial amount of storage. Therefore, we use a threefold folder structure:
1. The `damagesextremes` repo folder, which hosts i) all code, ii) small input files (e.g., country-level GDP projections from the SSP Database), and iii) intermediate data objects (e.g., RDS files). The latter, however, are not pushed to the repo due to their size
2. A `sharepoint` folder on the same hierarchy level as the 'damagesextremes' repo, which hosts i) larger input files (e.g., geospatial data on administrative boundaries from GADM) and ii) final outputs, such as figures and tables
3. External drives that host the ADM1-level aggregated climate indicator files and the GDP projections produced. Drives must be hard-mapped by the user in the scripts sourced by `src/01 Create Projections.R` and should provide sufficient storage (> 7.5TB) to host the resulting projection results

### Organization of the `damagesextremes` repository
* **data**: intermediate files created and saved out during the analysis and read in again, primarily to avoid expensive re-runs.
* **data/input**: all small input data
* **src**: all code used throughout the analysis
* **src/utils**: helper functions to load the climate data, produce GDP projections, calculate summary statistics of distributions efficiently, and create figures/tables. Functions in this folder are sourced in the header of the main scripts where required
* **src/projection_jobs**: stand-alone files that are sourced iteratively by `src/01 Create Projections.R` to create all necessary projections used in the paper

Note that the repo includes an `ccdamages-beyond-anntemp-v02.Rproj` file, which will automatically map the working directory to the `ccdamages-beyond-anntemp-v02` repo folder. If you do not use the Rproj file, you must set the working directory manually for all R scripts.

### Organization of the `sharepoint` directory
* **Data**: hosts all larger input files, such as administrative boundary shapefiles under `Data/GADM`
* **Figures**: contains all figures, which are created by the R scripts in the `ccdamages-beyond-anntemp-v02` repo
* **Tables**: contains all tables, which are created by the R scripts in the `ccdamages-beyond-anntemp-v02` repo

## System requirements
All scripts were executed on a high-performance cluster with 500GB RAM and 48 cores. Note that due to the large filesizes (particularly resulting from the Monte Carlo runs of the GDP projections), the scripts in their original version are likely to max out the memory of a conventional machine. Depending on the computational resources available to you, you might need to tweak the code to reduce memory requirements, e.g., via chunking parts of the analyses. External storage should be at least 7.5TB to host all GDP projection files for 1,000 Monte Carlo draws for damage function parameters.

## Scripts
To reproduce the results in the paper, run the following scripts in the specified order:
1. `src/create_gdp_weights_adm1.R`: uses gridded GDP data to calculate ADM1-level weights for aggregating data from ADM1 regions to the ADM0 level using GDP weighting. Runtime is about half a day on a high-performance cluster. Resulting file is saved out under `ccdamages-beyond-anntemp-v02/data`
2. `src/prepare_weights_and_shapefiles.R`: adds countries that have no ADM1-level subnational region in the GADM database but do feature in the SSP Database (Maldives, Aruba, Kiribati) to the ADM1-level administrative boundary shapefile and other auxiliary files (GDP & area weights) to ensure that they do feature in our projections. Resulting files are saved out under `ccdamages-beyond-anntemp-v02/data`
3. `src/calculate_stagg_areaweights.R`: calculates area weights to aggregate the downscaled CMIP6 data (0.25deg, ERA5 grid) to ADM1-level subnational regions in the GADM database using the `stagg` package and GADM shapefiles. Area weights are saved out under `ccdamages-beyond-anntemp-v02/data`
4. `src/aggregate_g025_to_adm1_stagg.R`: performs the area-weighted aggregation using the weights from the previous script for all bias-corrected climate indicators under consideration as well as for the raw (= not bias-corrected) monthly precipitation totals, which are used to inform where bias-corrected monthly precipitation deviation values are capped. Aggregated values for all climate indicators are saved out as feather files on the external drives
5. `src/calculate_monthly_precip_deviation.R`: calculates the monthly precipitation deviation, which Kotz et al., 2022 calculate only at the ADM1 level and not at the grid cell level. Results are saved as feather files on the external drives
6. `src/00 Estimate Regressions.R`: estimates the main specification from Kotz et al., 2022 (using their replication code) as well as the permutations of their regression model we use throughout our paper. Saves out regression models, coefficients, and variance-covariance matrices as R objects (RDS files) under `ccdamages-beyond-anntemp-v02/data` and produces the regression tables in the paper
7. `src/00 Prepare Analysis.R`: carries out multiple tasks to prepare the GDP projection calculation, such as producing overview tables for the Global Warming Level windows used or interpolating the SSP data on GDP. Resulting objects are saved out as RDS files under `ccdamages-beyond-anntemp-v02/data`
8. `src/01 Create Projections.R`: creates GDP projections at the ADM1 level for each CMIP6 model run for all years between the start of the +0.84C warming level period (our baseline) and 2100 for 1,000 Monte Carlo draws of dose-response function parameters, aggregates them up to the ADM0 level using GDP weighting, and saves them out as `.feather` files on the external drive. Results excl. Monte Carlo draws for dose-response function parameters are calculated for alternative specifications as well, such as including heat (Tx5d) or using an alternative way of bias-correcting CMIP6 outputs. Different steps involved in creating the projections are spread over different sub-scripts, which are located in `src/projection_jobs` and sourced by `src/01 Create Projections.R`. These sub-scripts also include calculating summary statistics for the GDP impact distributions and carrying out variance decomposition calculations. The combined runtime of these sub-scripts is multiple days even on a high-performance cluster
9. `src/03 Analyse Projections.R`: creates all figures, tables, and numeric results in the paper as well as the Supplementary Information except for our conceptual Figure 1 and the SI figure on dose-response functions. The script also includes analysis steps for variance decomposition and to calculate marginal effects, such that runtime is non-negligible as well. All figures and visualizations are saved out under the `sharepoint` directory
10. `src/04 Analyse Dose-Response.R`: creates the conceptual Figure 1 and the SI figure on dose-response functions, saving them out under the `sharepoint` directory

## Data files from external sources (timestamp = download date)
* `data/input/220803 SSP-Database GDP PPP.xlsx`: taken from the SSP Database
* `data/input/230108 World Bank Country Classification.xlsx`: taken from the World Bank Data Help Desk on the World Bank Country and Lending Groups
* `data/input/230503 World Bank GDP_2015USD.csv`: taken from the World Bank's World Development Indicator database
* `231030 World Bank WDI Total Population.xlsx`: population data taken from the World Bank's World Development Indicator database with additional CIA Factbook values added for Taiwan (see sheet ReadMe), which we use to calculate population exposures to tail risks
* `data/input/PT_master.csv`: the main dataset used by Kotz et al. (2022, Nature), which we use in `src/00 Estimate Regressions.R` to reproduce their regression model and estimate permutations
* `data/input/CallahanMankin2022_extremes_growth_panel_monthedd_1979-2016.csv`: the main dataset used by Callahan & Mankin (2022, Science Advances), from which we take their PDSI measure as well as their Tx5d values for regression model permutations
* `data/input/API_NY.GDP.DEFL.ZS_DS2_en_excel_v2_4353530.csv`: the US GDP deflator data used by Callahan & Mankin, which we deploy to show that inflation adjustments do not alter the results of the main regression model by Kotz et al. (2022)
* `data/input/df_gwl.csv`: global warming level windows calculated following Batibeniz, Hauser, and Seneviratne (2023, Earth Syst Dyn)
* `data/input/df_modelscen.csv`: a CSV file hosting the model, scenario and realization identifier of all 199 model runs used in our paper

## sessionInfo()
```
R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: openSUSE Leap 15.5

Matrix products: default
BLAS/LAPACK: /usr/local/OpenBLAS-0.3.21/lib/libopenblas_haswellp-r0.3.21.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8      
 [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] metR_0.13.0             PBSmapping_2.73.2       raster_3.6-11           sp_1.5-1                ncdf4.helpers_0.3-6     data.table_1.14.6       ncdf4_1.19             
 [8] stagg_0.0.0.9000        ggridges_0.5.4          dtplyr_1.2.2            ggcorrplot_0.1.4        feather_0.3.5           furrr_0.3.1             future_1.29.0          
[15] sf_1.0-9                ggpubr_0.5.0            ggpattern_1.0.1         PupillometryR_0.0.4     rlang_1.1.1             janitor_2.2.0           rnaturalearthdata_0.1.0
[22] rnaturalearth_0.1.0     xtable_1.8-4            readxl_1.4.1            gsubfn_0.7              proto_1.0.0             stargazer_5.2.3         clubSandwich_0.5.8     
[29] plm_2.6-3               forcats_0.5.2           stringr_1.5.0           dplyr_1.1.3             purrr_1.0.2             readr_2.1.3             tidyr_1.3.0            
[36] tibble_3.2.1            ggplot2_3.4.4           tidyverse_1.3.2        

loaded via a namespace (and not attached):
 [1] googledrive_2.0.0   colorspace_2.0-3    ggsignif_0.6.4      ellipsis_0.3.2      class_7.3-20        snakecase_0.11.1    fs_1.5.2            rstudioapi_0.15.0   proxy_0.4-27       
[10] listenv_0.8.0       fixest_0.11.1       fansi_1.0.3         lubridate_1.9.0     xml2_1.3.3          codetools_0.2-18    cachem_1.0.6        Formula_1.2-4       jsonlite_1.8.3     
[19] broom_1.0.5         dbplyr_2.2.1        compiler_4.2.2      httr_1.4.4          backports_1.4.1     fastmap_1.1.0       assertthat_0.2.1    gargle_1.2.1        cli_3.6.1          
[28] tools_4.2.2         gtable_0.3.1        glue_1.6.2          dreamerr_1.3.0      Rcpp_1.0.10         carData_3.0-5       cellranger_1.1.0    vctrs_0.6.4         nlme_3.1-160       
[37] lmtest_0.9-40       globals_0.16.2      rbibutils_2.2.13    rvest_1.0.3         timechange_0.1.1    lifecycle_1.0.3     collapse_1.8.9      rstatix_0.7.2       googlesheets4_1.0.1
[46] terra_1.7-39        MASS_7.3-58.1       zoo_1.8-11          scales_1.2.1        miscTools_0.6-28    hms_1.1.2           parallel_4.2.2      sandwich_3.0-2      memoise_2.0.1      
[55] bdsmatrix_1.3-6     stringi_1.7.8       checkmate_2.1.0     e1071_1.7-12        Rdpack_2.4          pkgconfig_2.0.3     lattice_0.20-45     tidyselect_1.2.0    parallelly_1.32.1  
[64] plyr_1.8.8          magrittr_2.0.3      R6_2.5.1            generics_0.1.3      DBI_1.1.3           pillar_1.9.0        haven_2.5.1         withr_2.5.0         units_0.8-0        
[73] abind_1.4-5         modelr_0.1.10       crayon_1.5.2        car_3.1-1           KernSmooth_2.23-20  utf8_1.2.2          tzdb_0.3.0          maxLik_1.5-2        grid_4.2.2         
[82] reprex_2.0.2        digest_0.6.30       classInt_0.4-8      numDeriv_2016.8-1.1 munsell_0.5.0       tcltk_4.2.2
```





