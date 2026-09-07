
rm(list=ls())

library('knitr')
library('tidyverse')
library('ggplot2')
library('scales')
library('dplyr')
library('png')
library('data.table')


source('stage2.R')
source('tmle.R')
source('aps.R')
source('incidence_calcs.R')

# UNBLINDED
dt <- readRDS("../endpoint_analytic_dataset_unblinded.rds")

recency_var <- 'recent_primary'
one_sided_test <- T

set.seed(1)
# primary analysis
inc_pri <- EndpointAnalysis(dt=dt, recency_var=recency_var,
                            save_stage_1=T)

# with secondary RITA parameters
inc_sec <- EndpointAnalysis(dt=dt, recency_var=recency_var,
                            parameter_pri=F)

# subgroups
inc_stable <- EndpointAnalysis(dt=dt[dt$stable==1,], recency_var=recency_var)

inc_men <- EndpointAnalysis(dt=dt[dt$male==1,], recency_var=recency_var)
inc_women <- EndpointAnalysis(dt=dt[dt$male==0,], recency_var=recency_var)

inc_young <- EndpointAnalysis(dt=dt[dt$young==1,], recency_var=recency_var)
inc_old <- EndpointAnalysis(dt=dt[dt$young==0,], recency_var=recency_var)

inc_kenya <- EndpointAnalysis(dt=dt[dt$kenya==1,], recency_var=recency_var)
inc_uganda <- EndpointAnalysis(dt=dt[dt$kenya==0,], recency_var=recency_var)

# post hoc subgroup - unhealthy alcohol use
inc_unhealthy <- EndpointAnalysis(dt=dt[which(alcohol_unhealthy==1),], recency_var=recency_var)

# incidence estimates when using TMLE to adjust for potentially differential ascertainment
library('ltmle')
# impute missing data on census variables
dt[is.na(dt$residence_months), 'residence_months'] <-  median(dt$residence_months, na.rm = T)
dt[is.na(dt$residence_nights), 'residence_nights'] <-  median(dt$residence_nights, na.rm = T)

# run TMLE
inc_tmle <- EndpointAnalysis(dt=dt, recency_var=recency_var, 
                             empirical_prev=F)

# incidence estimates when forced adjustment (post hoc)
inc_prev <- EndpointAnalysis(dt=dt, recency_var=recency_var,
                             verbose=T,
                             do.unadjusted=F, QAdj='prev_2yr')
inc_anc <- EndpointAnalysis(dt=dt, recency_var=recency_var, 
                            verbose=T,
                            do.unadjusted=F, QAdj='district_anc_prev')

# TMLE with Adaptive Pre-specification
inc_aps <- EndpointAnalysis(dt=dt, recency_var=recency_var, 
                            verbose=T,
                            do.unadjusted=F, QAdj='U')

inc_compile <- data.frame(rbind(overall=inc_pri, kenya=inc_kenya, 
                                uganda=inc_uganda, 
                                youth=inc_young, 
                                old=inc_old,
                                men=inc_men, 
                                women=inc_women))


inc_sens <- data.frame(rbind(stable=inc_stable, 
                             rita_sec=inc_sec, 
                             tmle_miss=inc_tmle, 
                             tmle_aps = inc_aps,
                             tmle_prev=inc_prev,
                             tmle_anc=inc_anc
                              ))
# save
save(inc_pri, inc_sec, inc_stable,
     inc_men, inc_women, inc_young, inc_old, inc_kenya, inc_uganda, inc_unhealthy,
     inc_tmle, inc_prev, inc_anc, inc_aps,
     inc_compile, inc_sens,
     file='incidence_unblinded.rdata')

