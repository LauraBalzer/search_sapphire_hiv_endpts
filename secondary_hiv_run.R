rm(list=ls())

library(data.table)
library(ltmle)
library(SuperLearner)

source('stage2.R')
source('tmle.R')
source("aps.R")
source('secondary_hiv_functions.R')

if (file.exists("stage1_secondary_hiv.csv")) {
  file.remove("stage1_secondary_hiv.csv")
}

set.seed(12345)

dt <- readRDS("../endpoint_analytic_dataset_unblinded.rds")
file_name <- 'secondary_unblinded.rdata'

DT_COMM <- unique(dt[, .(
  kenya,
  community_code,
  pair,
  studyarm,
  district_anc_prev
)])
stopifnot(DT_COMM[, .N] == 16)


# impute missing data on census variables
dt[is.na(dt$residence_months), 'residence_months'] <-  median(dt$residence_months, na.rm = T)
dt[is.na(dt$residence_nights), 'residence_nights'] <-  median(dt$residence_nights, na.rm = T)


dt <- dt[, .(
  subjid,
  community_code,
  studyarm,
  stable,
  age,
  male,
  residence_nights,
  residence_months,
  # needed for subgroups
  young, 
  kenya,
  # cascade
  hiv_prevalent, # identical to HIV
  hiv,
  pdx,
  art,
  vs,
  # viremia
  viremia_mp,
  viremia_10k, # this is really >99,999 ****
  # biomedical HIV prevention coverage 
  coverage,
  coverage_bin, # uptake 6mo
  prep_pep_use_2yrs, # uptake 2yrs
  coverage_2yr, # coverage 6mo
  # post-hoc subgroups
  pregnant, 
  alcohol,
  alcohol_unhealthy
  )
]


###############################################
## Run all endpoints
###############################################

one_sided <- T
code_qc <- F
verbose <- T

subgroups <- c('ALL', 'STABLE', 'MEN', 'WOMEN',
               'YOUNG', 'OLDER', 'KENYA', 'UGANDA',
               # following are post-hoc
               'NO_ALCOHOL', 'ALCOHOL_USE', 
               'NOT_UNHEALTHY', 'UNHEALTHY',
               'PREGNANT')

# First 90: P(DX=1 | HIV+)
pdx_hiv <- analyze_wrapper(dt=dt,
                           numerator_col  = "pdx",
                           denominator_col = "hiv",
                           one_sided=one_sided,
                           alt_smaller   = FALSE,  # H1: intervention > control
                           code_qc=code_qc, verbose=verbose, subgroups=subgroups
                           )

# Second 90: P(ART=1 | DX=1, HIV+)
art_pdx <- analyze_wrapper(dt=dt,
                           numerator_col  = "art",
                           denominator_col = "pdx",
                           one_sided=one_sided,
                           alt_smaller   = FALSE,   # H1: intervention > control
                           code_qc=code_qc, verbose=verbose, subgroups=subgroups
                           )


# Third 90: P(VS=1 | ART=1, DX=1, HIV+)
vs_art <- analyze_wrapper(dt=dt,
                          numerator_col  = "vs",
                          denominator_col = "art",
                          one_sided=one_sided,
                          alt_smaller   = FALSE,   # H1: intervention > control
                          code_qc=code_qc, verbose=verbose, subgroups=subgroups
                          
)

## Population-level VS among PWH 
vs_hiv <- analyze_wrapper(dt=dt,
                          numerator_col  = "vs",
                          denominator_col = "hiv",
                          one_sided=one_sided,
                          alt_smaller   = FALSE,   # H1: intervention > control
                          code_qc=code_qc, verbose=verbose, subgroups=subgroups
)


#*************************
#* VIREMIA 

# Population-level viremia (lower is better in intervention)
## Maya P’s simple identification result & approach
viremia <- analyze_wrapper(dt=dt,
                           numerator_col  = "viremia_mp",
                           one_sided=one_sided,
                           alt_smaller   = TRUE,    # H1: intervention < control
                           code_qc=code_qc, verbose=verbose, subgroups=subgroups
                           
)


# VL> 99,999
viremia_10k <- analyze_wrapper(dt=dt,
                               numerator_col  = "viremia_10k",
                               one_sided=one_sided,
                               alt_smaller   = TRUE,    # H1: intervention < control
                               code_qc=code_qc, verbose=verbose, subgroups=subgroups
                               
)


#*************************
#* PREVENTION COVERAGE 

# SUBSET TO TARGET POP: PERSONS NOT KNOWN TO HAVE HIV
dt_coverage <- dt[ dt$hiv==0 | is.na(dt$hiv),]
dim(dt_coverage)

## Prevention coverage (higher is better in intervention)
coverage <- analyze_wrapper(dt=dt_coverage,
                           numerator_col  = "coverage",
                           one_sided=one_sided,
                           alt_smaller   = FALSE,   # H1: intervention > control
                           code_qc=code_qc, verbose=verbose, subgroups=subgroups
                           )

## Note: the following were pre-specified 
## but not included in the original code upload
coverage_bin <- analyze_wrapper(dt=dt_coverage,
                                numerator_col  = "coverage_bin",
                                one_sided=one_sided,
                                alt_smaller   = FALSE,   # H1: intervention > control
                                code_qc=code_qc, verbose=verbose, subgroups=subgroups
)

prep_pep_use_2yrs <- analyze_wrapper(dt=dt_coverage,
                                numerator_col  = "prep_pep_use_2yrs",
                                one_sided=one_sided,
                                alt_smaller   = FALSE,   # H1: intervention > control
                                code_qc=code_qc, verbose=verbose, subgroups=subgroups)

coverage_2yr <- analyze_wrapper(dt=dt_coverage,
                            numerator_col  = "coverage_2yr",
                            one_sided=one_sided,
                            alt_smaller   = FALSE,   # H1: intervention > control
                            code_qc=code_qc, verbose=verbose, subgroups=subgroups)



# SAVE EVERYTHING

save(pdx_hiv, art_pdx, vs_art, vs_hiv,
     viremia, viremia_10k,
     coverage, coverage_bin,  prep_pep_use_2yrs, coverage_2yr,
     file=file_name)

