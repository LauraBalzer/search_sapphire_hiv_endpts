# code to run treatment-blind simulations to verify Type-I error control

rm(list=ls())
library('data.table')
source('stage2.R')
source('tmle.R')

source('incidence_calcs.R')

dt <- readRDS("../endpoint_analytic_dataset_unblinded.rds")

blind_me <- function(dt){ 
  # cat("Blinded analysis: assigning study arms at random\n")
  community_code_set <- dt[, unique(community_code)]
  # pick 8 communities at random to be intervention
  intervention_communities <- sample(community_code_set, 8)
  # print(sort(intervention_communities))
  dt[, studyarm := as.numeric(community_code %in% intervention_communities)]
  dt
}

J <-  5000

set.seed(12.16)

OUT <- data.frame(matrix(NA, nrow=J, ncol=4))
colnames(OUT) <- c('rel_pri', 'abs_pri', 'rel_sec', 'abs_sec')

for(j in 1:J){
  
  dt <- blind_me(dt)
  # dt[, .N, keyby = .(studyarm, community_code)]

  # primary analysis
  temp <- EndpointAnalysis(dt=dt, recency_var='recent_primary',
                           parameter_pri=T,
                           one_sided_test=T)
  OUT[j, 'rel_pri'] <- temp[,'Rel.pval']
  OUT[j, 'abs_pri'] <- temp[,'Abs.pval']
  
  # secondary analysis with secondary RITA parameters
  temp <- EndpointAnalysis(dt=dt, recency_var='recent_primary',
                           parameter_pri=F,
                           one_sided_test=T)
  OUT[j, 'rel_sec'] <- temp[,'Rel.pval']
  OUT[j, 'abs_sec'] <- temp[,'Abs.pval']
  
  print(j)
}

colMeans(OUT < 0.05)


save(OUT,
        file=paste0('inc_txt_blind_reps', J,  '.rdata') )
