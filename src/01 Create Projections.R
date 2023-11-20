rm(list = ls())

# NOTE: since creating projections is very computation-intensive, each part is
# outsourced into a stand-alone script under src/projection_jobs. The script here
# serves as a master file to clarify the order of the jobs. Depending on the
# specs of the machine used, these jobs might max out the memory or take a very
# long time to compute. Total runtime on a machine with 500GB RAM and 48 workers
# was a number of days and will increase substantially if a lower number of
# workers is available.

# ADM1-, ADM0- and global projections using point estimates for dose-response functions:
source(file.path("src", "projection_jobs", "create_df_adm1_bc.R"))

# ADM1-, ADM0- and global projections using point estimates for dose-response functions when including Tx5d
source(file.path("src", "projection_jobs", "create_df_adm1_bc_tx5d.R"))

# ADM1-, ADM0- and global projections using point estimates for dose-response functions for single model run w/ alternative bias correction
source(file.path("src", "projection_jobs", "create_df_adm1_bc_altbiascor.R"))

# ADM0-level projections using full Monte Carlo distribution for dose-response functions (= main results)
source(file.path("src", "projection_jobs", "create_df_adm0_bc_fulldistr.R"))

# compiling all ADM0-level projections using full MC distribution for each sovereign country
# NOTE: this is required because due to memory constraints, we cannot simply read in all ADM0 results at once
source(file.path("src", "projection_jobs", "compile_adm0_bc_fulldistr.R"))

# aggregate all ADM0-level projection using full MC distribution for dose-response functions (= main results) to global level
source(file.path("src", "projection_jobs", "aggregate_all_fulldistr_to_global.R"))

# calculate summary statistics for the global distribution with and without Monte Carlo draws for dose-response functions
# NOTE: this also calculates weights for each CMIP6 model used to ensure that all models have same sampling probability
source(file.path("src", "projection_jobs", "summarise_global_bc_fulldistr.R"))

# calculate summary statistics at ADM0 level for distribution with Monte Carlo draws for dose-response functions
source(file.path("src", "projection_jobs", "summarise_adm0_bc_fulldistr.R"))

# perform variance decomposition following Hsiang et al. 2017 for global distribution with Monte Carlo draws for dose-response functions
source(file.path("src", "projection_jobs", "decompvar_global_bc_fulldistr.R"))
