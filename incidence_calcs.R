# main functions to analyze the primary endpoint: HIV incidence
# code written by Josh Schwab
# documented & edited by Laura Balzer

EndpointAnalysis <- function(dt, 
                             recency_var='recent_primary', 
                             empirical_prev=T,
                             parameter_pri=T, 
                             one_sided_test=T, 
                             verbose=F,
                             save_stage_1=F, 
                             force_adj=F, QAdj=NULL) {
  
  
  stopifnot(recency_var %in% c("recent_primary", "recent_stats", 
                               "recent_no_vl", "recent_all_ages"))
  
  # get community-specific counts of 
  dt_analysis <- dt[, .(.N, 
                        n_hivstatus = sum(!is.na(hiv_prevalent)),
                        n_hivpos = sum(hiv_prevalent == 1, na.rm = T),
                        n_elig = sum(eligible_stats == 1),
                        n_recent = sum(get(recency_var) == 1, na.rm = T),
                        # not needed for incidence calcs, but proxy for prev at baseline
                        n_dx_2yr = sum(hiv_dx_old==1, na.rm=T)), 
                    by = .(community_code, studyarm, kenya, pair, district_anc_prev)]
  
  
  # PrevH = Prevalence of HIV 
  if(empirical_prev){
    # primary: unadjusted prevalence estimate
    # denominator is number with known status (not N)
    dt_analysis[, PrevH := n_hivpos / n_hivstatus] 
  } else{
    # sensitivity: TMLE to adjust 
    dt_prev <- dt[, .( PrevH = est_prev_tmle(.SD)), by = community_code]
    
    dt_analysis <- merge(
      dt_analysis,
      dt_prev,
      by = "community_code"
    )
  }


  
  # PrevR = % 'recent' by biomarker assay among those with prevalent HIV
  dt_analysis[, PrevR := n_recent / n_hivpos]
  # with >99% coverage among PWH, skipping TMLE in sensitivity
  
  # print( dt_analysis[,c('N', 'n_hivstatus', 'n_hivpos', 'n_elig', 'n_recent', 'PrevH', 'PrevR')])
  
  # RITA parameters
  if(parameter_pri){
    this_mdri <- 130
    this_frr <- 0
  } else{
    this_mdri <- 198.90
    this_frr <- 0.0155
  }
  
  dt_analysis[, MDRI := this_mdri] # 
  # community specific scaling of the FRR
  # Not eligible if deterministically not recent (>2yr, on ART, VL<75)
  # P(FR) = P(FR | eligible)*P(eligible) + P(FR | not-eligible)*P(not-eligible)
  # Second term ==0 
  dt_analysis[, FRR := this_frr * n_elig / n_hivpos]  
  dt_analysis[, BigT := 730] #constant
  
  # Per-community incidence
  dt_analysis[, inc_hat := EstIncidence(PrevH, PrevR, MDRI, FRR, BigT), 
              by = community_code]
  
  # ensure that 0 recent = 0 incidence
  dt_analysis[dt_analysis$n_recent==0, 'inc_hat'] <- 0
  
  if(verbose){
    print(paste0('Communities analyzed:', nrow(dt_analysis)))
    print(paste0('With zero recent:', sum(dt_analysis$inc_hat==0)))
  }
  
  # OUTPUT STAGE 1 
  if(save_stage_1){
    dt_out <- dt_analysis[,c('kenya', 'pair','community_code','studyarm',
                             'inc_hat')]
    write.csv(dt_out,'stage1_hiv_inc.csv', row.names = F)
  }
  
  # STAGE 2 
  
  # primary analysis is a one-sided test
  if(one_sided_test){
    one.sided <- T
    alt.smaller <- T
  } else{
    one.sided <- F
    alt.smaller <- NA
  }
  
  arm_cols <-  c(paste0('Txt', c('.est', '.CI.lo', '.CI.hi')),
                 paste0('Con', c('.est', '.CI.lo', '.CI.hi')) )
  
  effect_cols <-  c('est', 'CI.lo', 'CI.hi', 'pval')
  
  these_cols <- c(arm_cols, effect_cols)
  
  if(!force_adj){
    # primary: unadjusted effectiveness estimate (Stage 2)
    df <- as.data.frame(dt_analysis[, .(id = 1:16, U = 1, alpha = 1, 
                                        A = studyarm, Y = inc_hat)])
    do.unadjusted <- T
    QAdj <- NULL
  }else {
    
    if(QAdj=='prev_2yr'){
      # post hoc analysis adjusting for proxy of baseline prev
      dt_analysis[, prev_2yr := n_dx_2yr / N] 
      df <- as.data.frame(dt_analysis[, .(id = 1:16, U = 1, alpha = 1, 
                                          W = prev_2yr,
                                          A = studyarm, Y = inc_hat)])
      # write.csv(df,'stage1_hiv_inc_prev.csv', row.names = F)
      
    } else{
      # post hoc analysis adjusting for ANC prev at baseline
      df <- as.data.frame(dt_analysis[, .(id = 1:16, U = 1, alpha = 1, 
                                          W = district_anc_prev,
                                          A = studyarm, Y = inc_hat)])
      # write.csv(df,'stage1_hiv_inc_anc.csv', row.names = F)
      
    }
   
    print(paste(QAdj, '; int: ', round(mean(df[df$A==1,'W'])*100,1), '%; ',
                'soc: ', round(mean(df[df$A==0,'W'])*100,1), '%.'
    ))
    
    do.unadjusted <- F
    QAdj <- 'W'
    verbose <- T
    
  }
  
  
  # for the relative effect
  rel <- Stage2(goal='aRR', data.input=df, do.unadjusted=do.unadjusted, 
                QAdj=QAdj,
                one.sided=one.sided, alt.smaller=alt.smaller,
                verbose=verbose)
  rel <- rel[, these_cols]
  
  # for the absolute effect
  abs <- Stage2(goal='RD', data.input=df, do.unadjusted=do.unadjusted, 
                QAdj=QAdj,
                one.sided=one.sided, alt.smaller=alt.smaller,
                verbose=verbose)
  abs <- abs[, these_cols]
  
  yay <- cbind(abs, rel[,effect_cols])
  colnames(yay) <- c(arm_cols, paste0('Abs.', effect_cols), paste0('Rel.', effect_cols))
  
  yay
}



EstIncidence <- function(PrevH, PrevR, MDRI, FRR, BigT) {
  # code adapted from https://github.com/SACEMA/inctools/blob/master/R/incidence.R
  # with relevant line numbers below
  # See `compare_incidence_equations.R` for alternative implementations
  # All give the same answer
  
  # SIMPLE FORMULA FOR INCIDENCE - line 79 
  I_EST <- function(prevH, prevR, mdri, frr, bigt) {
    prevH * (prevR - frr)/((1 - prevH) * (mdri - frr * bigt))
  }
  
  # NOTE: mdri and bigt must be in years - line 297
  # Convert to years
  MDRI <- MDRI/365.25
  BigT <- BigT/365.25
  I_Est <- I_EST(prevH = PrevH, prevR = PrevR, mdri = MDRI, frr = FRR, bigt = BigT)
    
  return(I_Est)
}



community_counts <- function(dt) {
  dt_analysis <- dt[, .(
    .N,
    n_participant = sum(!census_only),
    n_hivstatus = sum(!is.na(hiv_prevalent)),
    n_hivpos = sum(hiv_prevalent == 1, na.rm = TRUE),
    n_elig = sum(eligible_stats == 1),
    n_recent = sum(recent_primary == 1, na.rm = TRUE)
  ), by = .(community_code, studyarm)]
  
  # Sum across communities within each study arm
  dt_arm <- dt_analysis[, lapply(.SD, sum), 
                        .SDcols = c("N","n_participant", "n_hivstatus", "n_hivpos", "n_elig", "n_recent"), 
                        by = studyarm]
  
  list(by_community = dt_analysis, by_studyarm = dt_arm)
}

est_prev_tmle <- function(dt_comm){
  
  # adjust for age, sex, mobility, 
  dt_in <- dt_comm[, .(male, age, residence_nights, residence_months,
                       C = BinaryToCensoring(is.censored = is.na(get('hiv'))), 
                       Y = get('hiv'))]
  suppressMessages({
    fit <- ltmle(
      data   = dt_in,
      Anodes = NULL, 
      abar = NULL, 
      Cnodes = "C",
      Ynodes = "Y",
      SL.library = 'glm',  
      variance.method = "ic",
      estimate.time = FALSE
    )
  })
  unname(fit$estimates["tmle"])
}

