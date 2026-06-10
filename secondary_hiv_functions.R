# main functions to analyze the secondary HIV-related endpoints
# code written by Josh Schwab
# documented & edited by Laura Balzer

analyze_wrapper <- function(dt, numerator_col, denominator_col = "none",
                           one_sided=T, alt_smaller, code_qc=F, verbose=F) {
 
  OUT <- NULL 

  if(code_qc){
    subgroups <-  'ALL'
  } else{
    subgroups <- c('ALL', 'STABLE', 'MEN', 'WOMEN',
                   'YOUNG', 'OLDER', 'KENYA', 'UGANDA')
  }
 
  for(j in subgroups){
    
    dt_in <- get_subgroups(dt, j)
    if(!code_qc){ print(paste0(j, '- participants: ', nrow(dt_in))) }
    
    append_stage1_csv <- j == "ALL" & numerator_col %in% 
      c("coverage", "pdx", "art", "vs", "viremia_mp", "viremia_10k")   
    
    # STAGE 1: TMLE
    # STAGE 2: TMLE & UNADJ
    temp <- analyze(dt_in=dt_in,
                   numerator_col = numerator_col,
                   denominator_col = denominator_col,
                   unadj_stage1=FALSE,
                   one_sided=one_sided,
                   alt_smaller = alt_smaller, 
                   verbose=(verbose & j=='all'), 
                   code_qc=code_qc,
                   append_stage1_csv = append_stage1_csv
                   )
    rownames(temp) <- paste0(j,'_t_', rownames(temp))
    OUT <- rbind(OUT, temp)  
    
    # STAGE1: UNADJ
    # STAGE2: TMLE & UNADJ
    temp <- analyze(dt_in=dt_in,
                    numerator_col = numerator_col,
                    denominator_col = denominator_col,
                    unadj_stage1=TRUE,
                    one_sided=one_sided,
                    alt_smaller = alt_smaller, 
                    verbose=(verbose & j=='all'), 
                    code_qc=code_qc,
                    append_stage1_csv = F)['u',]
    rownames(temp) <- paste0(j,'_u_', rownames(temp))
    OUT <- rbind(OUT, temp)  
    
  }
  
  OUT
  
}




# analyze- run Two-Stage TMLE
# Stage 1: estimate the endpoint for each community
#   if a non-ratio estimand, then runs TMLE once for each community
#   if a ratio estimand, then runs TMLE twice for each community
#   note: .SD subsets dt on community 
# Stage 2: evaluate effectiveness
#   apply TMLE with APS and unadjusted estimator
#   incorporates one-sided hypo tests

analyze <- function(dt_in, numerator_col, denominator_col = "none", unadj_stage1=F,
                    one_sided=T, alt_smaller, verbose=F, code_qc=F, append_stage1_csv) {
  
  dt <- copy(dt_in)
  
  if (denominator_col == "none") {
    # for all non-ratio estimands (e.g. coverage)
    comm_est <- dt[, .(est = run_tmle(.SD, Ynode = numerator_col, unadj_stage1=unadj_stage1)), 
                   keyby = community_code]
  } else {
    # for all ratio-type estimands (i.e. HIV care cascade)
    comm_est <- dt[, .(est = run_tmle(.SD, Ynode = numerator_col, unadj_stage1=unadj_stage1) / 
                         run_tmle(.SD, Ynode = denominator_col, unadj_stage1=unadj_stage1)), 
                   keyby = community_code]
  }
  
  # subset only on community code and community's endpt estimate
  comm_est <- comm_est[, .(community_code, est)]

  # if(!code_qc){ print(format_round(comm_est$est, 3)) }

  # merge in cluster-level variables
  dt_comm <- merge(
    DT_COMM,
    comm_est,
    by = "community_code"
  )

  stopifnot(dt_comm[, .N] %in% c(8,16) ) 
  # if(!code_qc){ print(paste0('# communities: ', dt_comm[, .N])) }
  
  #error if any NA values in any column of dt_comm
  stopifnot(!any(is.na(dt_comm)))
  
  # build Stage2 input with covariates
  df <- as.data.frame(dt_comm[, .(id = 1:.N, 
                                  alpha = 1,
                                  U = 1,
                                  district_anc_prev,
                                  A = studyarm,
                                  Y = est
  )])
  if(!code_qc){print(mean(df$Y))}
  
  if (append_stage1_csv) {
    dt_csv <- dt_comm[, .(numerator_col, denominator_col, kenya,
                           pair, community_code, studyarm, est)]
    stopifnot(!any(is.na(dt_csv)))
    fwrite(dt_csv, file = "stage1_secondary_hiv.csv", append = T)
  }

  # simple comparison of mean outcomes
  unadj <- run_stage2(df=df,
                      unadj_stage2=T, 
                      one_sided=one_sided,
                      alt_smaller = alt_smaller,
                      verbose=F) 
  
  if(unadj_stage1){
    tmle <- rep(NA, length(unadj))
  }else{
    # TMLE to improve efficiency
    tmle <- run_stage2(df=df,
                       unadj_stage2=F, 
                       one_sided=one_sided,
                       alt_smaller = alt_smaller,
                       verbose=verbose)
  }
  
  
  data.frame(rbind(t=tmle, 
                   u=unadj))
}


# run_tmle: calls ltmle to estimate community endpts in Stage 1
#   adjust for differences between those measured/missed
#   these covariates are limited to data collected via census (age, sex, mobility)

run_tmle <- function(data_in, Ynode, unadj_stage1=FALSE) {
  
  # most estimands do not have deterministic knowledge
  det_Q_fun <- NULL
  
  # subset on the relevant columns
  # generate a censoring node to indicate missingness on the outcome
  # note: get() grabs column corresponding to Ynode
  
  if(Ynode %in% c('pdx', 'art') | unadj_stage1){
    # No missingness - no adjustment for prior dx & ART use
    # sensitivity analyses with empirical proportion (unadj_stage1=T)
    data <- data_in[, .(C = BinaryToCensoring(is.censored = is.na(get(Ynode))), 
                        Y = get(Ynode))]
  } else{
    # adjust for age, sex, mobility, 
    data <- data_in[, .(male, age, residence_nights, residence_months,
                        C = BinaryToCensoring(is.censored = is.na(get(Ynode))), 
                        Y = get(Ynode))]
  }

  
  if (Ynode == "vs") {
    # P(vs^*=1,art=1,pdx=1,HIV^*=1) = sum_w P(vs=1|Δ_vl=1,art=1,w)×P(W=w|art=1)×P(art=1)
    # use deterministic knowledge that must be on ART if have VS
    data <- cbind(data_in[, .(art)], data)
    det_Q_fun <- function(data_det, current.node, nodes, called.from.estimate.g) {
      is.deterministic <- data_det$art == 0
      Q.value <- 0
      return(list(is.deterministic = is.deterministic, Q.value = Q.value))
    }
  } 
  
  suppressMessages({
    fit <- ltmle(
      data   = data,
      Anodes = NULL, 
      abar = NULL, 
      Cnodes = "C",
      Ynodes = "Y",
      SL.library = 'glm',  
      deterministic.Q.function = det_Q_fun,
      variance.method = "ic",
      estimate.time = FALSE
    )
  })
  # return TMLE estimate
  unname(fit$estimates["tmle"])
}




run_stage2 <- function(df,
                       unadj_stage2=F,
                       one_sided,
                       alt_smaller, 
                       verbose) {
  
  if(nrow(df)==16 & (mean(df$Y)>0.02) ){
    # specify candidate adjustment set
    # need to include the unadjusted 
    cand_adj <- c("U", "district_anc_prev")
  } else{
    # using unadjusted for country subgroups and rare outcomes
    unadj_stage2 <- T
  }
  
  # absolute/additive effect
  abs <- Stage2(goal = 'RD', target='clust',
                data.input   = df,
                do.unadjusted = unadj_stage2,
                one.sided    = one_sided,
                alt.smaller  = alt_smaller,
                do.data.adapt = !unadj_stage2,
                cand.QAdj = cand_adj,
                gAdj= 'U', # pre-specify unadjusted pscore
                verbose=verbose
  )
  
  # relative effect
  rel <- Stage2(goal = 'aRR', target='clust',
                data.input   = df,
                do.unadjusted = unadj_stage2,
                one.sided    = one_sided,
                alt.smaller  = alt_smaller,
                do.data.adapt = !unadj_stage2,
                cand.QAdj = cand_adj,
                gAdj= 'U', # pre-specify unadjusted pscore
  )

  arm_cols   <-  c(paste0('Txt', c('.est', '.CI.lo', '.CI.hi')),
                   paste0('Con', c('.est', '.CI.lo', '.CI.hi')) )
  effect_cols <-  c('est', 'CI.lo', 'CI.hi', 'se', 'pval', 'QAdj')
  these_cols  <- c(arm_cols, effect_cols)
  
  yay <- cbind(abs[,these_cols], rel[,effect_cols])
  colnames(yay) <- c(arm_cols, paste0('Abs.', effect_cols), 
                     paste0('Rel.', effect_cols))
  
  yay
}


get_subgroups <- function(dt, this_subgroup){
  
  if(this_subgroup=='ALL'){
    dt_in <- dt
  } else if(this_subgroup=='STABLE'){
    dt_in <- dt[stable == 1,]
  } else if(this_subgroup=='MEN'){
    dt_in <- dt[male == 1,]
  } else if(this_subgroup=='WOMEN'){
    dt_in <- dt[male == 0,]
  } else if(this_subgroup=='YOUNG'){
    dt_in <- dt[young == 1,]
  } else if(this_subgroup=='OLDER'){
    dt_in <- dt[young == 0,]
  } else if(this_subgroup=='KENYA'){
    dt_in <- dt[kenya == 1,]
  } else if(this_subgroup=='UGANDA'){
    dt_in <- dt[kenya == 0,]
  } 
  
  dt_in 
}
