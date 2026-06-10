##############
# stage2.R 
# R code to implement all Stage2 analyses to compare intervention effects between arms
# Using TMLE with and without Adaptive Pre-specification 

# Key citations on TMLE with Adaptive Pre-specification
#   Balzer at al. Stat Med. 2026: https://pubmed.ncbi.nlm.nih.gov/27436797/
#   Balzer et al. Biometrics. 2023: https://pubmed.ncbi.nlm.nih.gov/38446441/
#   Balzer et al. Clin Trials. 2026: https://pubmed.ncbi.nlm.nih.gov/41762058/

# Modified from Stage2_Functions.R in https://github.com/LauraBalzer/TwoStageTMLE
#   Paper: https://pubmed.ncbi.nlm.nih.gov/34939083/

# With edits from Josh Nugent (https://github.com/joshua-nugent):
#   to code 1-sided p-values & rescale outcomes to be bounded in [0,1]

# With lots of edits to make code more generalizable/robust
#---------------------------------------

# Stage2: Main function for estimation and inference 
# input: 
#	  goal: aRR= arithmetic risk ratio;  RD=risk difference; OR= odds ratio (not recommended)
#   target: cluster-level effect ("clust")  or pooled-indv effect ("indv")
#	  data.input: observed data 
#  
# --- options for covariate adjustment to improve precision
#   (1) unadjusted effect estimator: contrast in arm-specific means
#     do.unadjusted: if TRUE, overides other input (default=F)
#  
#   (2) prespecified but fixed adjustment 
#     QAdj: adjustment variables for the conditional mean outcome 
#     Qform: form of the working outcome regression 
#     gAdj: adjustment variables for the propensity score 
#     gform: form of the working pscore regression
#
#   (3) TMLE with Adaptive Pre-specification (APS) - RECOMMENDED
#     do.data.adapt: indicator to do APS
#     *specifying the candidate estimation approaches: see aps/get.cand.adj 
#     cand.QAdj: candidate adjustment variables for conditional mean outcome 
#     cand.Qform: candidate adjustment approaches for conditional mean outcome
#     cand.gAdj: candidate adjustment variables for the propensity score
#     cand.gform: candidate adjustment approaches for the propensity score
#     remove.pscore: indicator if a candidate adjustment variable should be removed for estimation of the pscore 
#           IF it was selected during estimation of the conditional mean outcome 
#           - recommended for trials with very few indpt units
#     * specifying the cross-validation (CV) scheme; defaults to LOOCV if # indpt units <=40
#     V: number of folds for cross-validation 
#
# --- 
# 
#   do.cv.variance: indicator to get cross-validated variance (only for effect estimates)
#   break.match: indicator to break matched pairs if pair-matched trial
#   one.sided: indicator of one-sided hypothesis test
#   alt.smaller: must specify the direction of the alternative hypothesis if one-sided test
#	  verbose: indicator to print updates 
#   psi:true value of effect, if known 
#   return.IC: indicator of whether to return the influence curve
#
# output: point estimate and inference

#-------------------
# This version 
# (1) transforms outcomes that are outside of [0,1] as in Chpt7 of TLB
# (2) one-sided hypothesis testing (requires specifying the direction of the alternatives)
# (3) generalizes cross-validation to V-fold if # indpt units > 40 
# (4) can adjust for multiple covariates. candidate must be a list
# (5) can estimate the outcome regression and pscore with more things than glm
#   current options: stepwise regression (with/without interactions), LASSO, 
#   multivariate adaptive regression splines (mars)
# (6) for CRTs, incorporates weights to target indv or cluster-level effect 
#   (regardless of level of data)
#   requires user-specified weight calculation (alpha)


#-------------------
# REQUIRES
# dummy column U=1

# pair: if pair-matched and want to keep pairs, the column indicating pairs

# column for weights (alpha)
#   set =1 for individually randomized trials
#   BUT for cluster randomized trials: 
#     value of the weights depends on the target of inference and data level 
#     Details in Benitez et al. Stat Med 2023
#     let J=number of clusters, N_j = cluster-specific sample size, N_tot = total # participants= sum_j N_j
#       if target='clust' with cluster-level data, alpha=1 
#       if target='clust' with indv-level data, alpha= 1/N_j
#       if target='indv' with cluster-level data, alpha= J/N_tot*N_j
#       if target='indv' with indv-level data, then alpha=1
#     for demonstration, see sim2.R in https://github.com/LauraBalzer/Comparing_CRT_Methods
#   weights should sum to the total # of randomized units
#   future work: Make this more general 
#   
#-------------------

Stage2 <- function(goal='aRR', target='clust', data.input, do.unadjusted=F, 
                   QAdj=NULL, Qform='glm', gAdj=NULL, gform='glm',
                   do.data.adapt =F, 
                   cand.QAdj=NULL, cand.Qform='glm', cand.gAdj=NULL, cand.gform='glm',
                   V=5, remove.pscore=F, do.cv.variance=F,
                   break.match=T, one.sided=F, alt.smaller=NULL, verbose=F, psi=NA,
                   return.IC=F){	
  
  #=====================================================
  # update: TRANSFORM the outcome as in Chpt7 of TLB 
  # no impact on outcomes already bounded in [0,1]
  if(max(data.input[,'Y']) > 1){
    scale_value <- max(data.input[,'Y'])
    # print(paste0('max Y: ', scale_value))
  } else {
    scale_value <- 1
  }
  if(min(data.input[,'Y']) < 0){
    scale_value_min <- min(data.input[,'Y'])
   #  print(paste0('min Y: ', scale_value))
  } else {
    scale_value_min <- 0
  }
  data.input[,'Y'] <- (data.input[,'Y'] - scale_value_min) / (scale_value - scale_value_min)

  #=====================================================
  # ADAPTIVE PRESPECIFICATION
  # update: flexibility in CV-scheme and candidate prediction algorithms
  if(do.data.adapt){
    # print('running APS')
    select <- do.adaptive.prespec(goal=goal, target=target, break.match = break.match, 
                                  Ldata= data.input, V=V,
                                  cand.QAdj=cand.QAdj, cand.Qform=cand.Qform,
                                  cand.gAdj=cand.gAdj, cand.gform=cand.gform,
                                  remove.pscore=remove.pscore,
                                  QAdj=QAdj, gAdj=gAdj,
                                  scale_value = scale_value, scale_value_min = scale_value_min,
                                  verbose=verbose)
    
    Q.index <- select$Q.index
    QAdj <- select$QAdj
    Qform <- select$Qform
    g.index <- select$g.index
    gAdj <- select$gAdj	
    gform <- select$gform
    
  } else{
    # QAdj <- gAdj <- 'U'
    Q.index <- g.index <- 1
  }
  
  # RUN FULL TMLE WITH ADJUSTMENT SET 
  if( do.unadjusted){
    # warning this does not incorporate weights
    # do not use if target=indv & data are at cluster-level
    est <- do_unadjusted(goal=goal, train=data.input,
                         scale_value = scale_value, scale_value_min = scale_value_min,
                         verbose=verbose)
  } else{
    # run TMLE
    est <- do.TMLE(goal=goal, target=target, train=data.input, QAdj=QAdj, Qform=Qform, 
                   gAdj=gAdj, gform=gform, 
                   scale_value = scale_value, scale_value_min = scale_value_min,
                   doing.CV=F, verbose=verbose)  
  }
  
  
  # GET INFERENCE 
  n.clust <- length(unique(data.input$id)) 
  if(verbose){ print(paste0('# indpt units: ', n.clust))}
  
  # Get point estimates of the treatment-specific mean
  R1 <- est$R1
  R0 <- est$R0
  
  # Note: this only gives standard (not cross-validated) inference
  Txt <- get.inference(psi.hat=R1, se=sqrt(est$var.R1), df=(n.clust-2))[,c('est','CI.lo','CI.hi','se')]
  Con <- get.inference(psi.hat=R0, se=sqrt(est$var.R0), df=(n.clust-2))[,c('est','CI.lo','CI.hi','se')]
  
  # Now: for the intervention effect 
  #  the point estimate on the relevant scale for getting inference
  if( goal=='aRR' ){
    psi.hat <- log(R1/R0)
  } else if (goal=='RD'){
    psi.hat <- R1- R0
  } else if (goal=='OR'){
    psi.hat <- log( R1/(1-R1)*(1-R0)/R0)
  }
  
  if(break.match){
    # if breaking the match, set df to (#clusters -2)
    df <- n.clust - 2
    var.hat <- est$var.break
  } else{
    # if preserving the match, set df to (#pairs-1)
    df <- length(unique(data.input$pair)) -1 
    var.hat <- est$var.pair
  }
  
  if(verbose){ print(paste0('hypo test is two-sided: ', !one.sided))}
  
  inference <- get.inference(goal=goal, psi=psi, psi.hat=psi.hat, se=sqrt(var.hat), df=df,
                             one.sided=one.sided, alt.smaller = alt.smaller)

  if(do.cv.variance){
    # if getting cross-validated inference
    inference.CV <- get.inference(goal=goal, psi=psi, psi.hat=psi.hat, se=sqrt(select$var.CV), df=df,
                                  one.sided=one.sided, alt.smaller = alt.smaller)
    
    est.df<-  data.frame(Txt=Txt, Con=Con, psi=psi, inference, CV=inference.CV, 
                         QAdj=Q.index, Qform=est$Qform, 
                         gAdj=g.index, gform=est$gform)
  } else{
    est.df <-  data.frame(Txt=Txt, Con=Con, psi=psi, inference, 
                          QAdj=Q.index, Qform=est$Qform, 
                          gAdj=g.index, gform=est$gform)
  }
  

  if(return.IC){
    RETURN <- list(IC=est, est.df=est.df)
  } else{
    RETURN <- est.df
  }
  RETURN
}

#-----------------------------------------------------#-----------------------------------------------------
# get.IC.variance - function to do influence curve-based variance estimate 


# input: 
#		goal: aRR= arithmetic risk ratio; RD for the risk difference; OR for the odds ratio
#   target: cluster-level effect ("clust") or pooled-indv effect ("indv")
#		Vdata: dataset 
#   R1: mean outcome under intervention
#   R0: mean outcome under SOC
#   sample.effect: indicator of targeting the sample (vs. population) estimand
#   scale_value: maximum value for outcome scaling
#   scale_value_min: minimum value for outcome scaling
#   doing.CV: indicator of getting cross-validation 
#   verbose: indicator ot print updates/output 
#
# output: estimated IC & variance - preserving/breaking the match
#   on log scale for if goal='aRR' or 'OR'

#-----------------------------------------------------#-----------------------------------------------------
get.IC.variance <- function(goal, target, Vdata, R1=NA, R0=NA, sample.effect=T,  
                            scale_value = 1, scale_value_min = 0, 
                            doing.CV=F, verbose=F){
  
  # number of randomized units
  J <- length(unique(Vdata$id))
  
  # calculate the relevant components of the IC 
  if(sample.effect){
    # default - assume interest is in the sample effect
    DY1 <- Vdata$alpha*Vdata$H.1W*(Vdata$Y - Vdata$Qbar1W.star)
    DY0 <- Vdata$alpha*Vdata$H.0W*(Vdata$Y - Vdata$Qbar0W.star)
  } else{
    # calculate the IC for population effect (extra term for DW)
    DY1 <- Vdata$alpha*( Vdata$H.1W*(Vdata$Y - Vdata$Qbar1W.star) + Vdata$Qbar1W.star - R1 )
    DY0 <- Vdata$alpha*( Vdata$H.0W*(Vdata$Y - Vdata$Qbar0W.star) + Vdata$Qbar0W.star - R0 )	
  }
  
  # unscale 
  DY1 <- DY1*(scale_value - scale_value_min) + scale_value_min
  DY0 <- DY0*(scale_value - scale_value_min) + scale_value_min
  
  # if individual-level data, then need to aggregate the IC to the cluster-level 
  # approach for aggregation depends on the target effect
  if( length(DY1) > J ) {
    
    if(target=='clust'){
      # Data are indv-level; target is cluster-level 
      if(!doing.CV) print('data=indv; target=clust')
      DY1 <- aggregate(DY1, by=list(Vdata$id), sum)[,-1]
      DY0 <- aggregate(DY0, by=list(Vdata$id), sum)[,-1]
      
    }else if(target=='indv') {
      # Data are indv-level; target is indv-level 
      if(!doing.CV & verbose) print('data=indv; target=indv')
      DY1 <- c(ltmle:::HouseholdIC(as.matrix(DY1), id = Vdata$id))
      DY0 <- c(ltmle:::HouseholdIC(as.matrix(DY0), id = Vdata$id))
      
    } 
    
    # for the pair-matched IC also need to aggregate to the cluster-level
    # Vdata <- aggregate(Vdata, by=list(Vdata$id), mean)[,-1]
  } 
  
  # INFLUENCE CURVES ARE NOW AT THE LEVEL OF THE RANDOMIZED UNIT
  if(goal=='RD'){
    # going after RD, easy IC
    DY <-  DY1 - DY0
    
  } else if (goal=='aRR'){ 
    # going after aRR, then get IC estimate on log scale
    #	i.e. Delta method for log(aRR) = log(R1) - log(R0)
    DY <- 1/R1*DY1 - 1/R0*DY0
    
  } else if(goal=='OR'){
    # Delta method for log(OR)
    DY <- 1/R1*DY1 + 1/(1-R1)*DY1 - 1/(1-R0)*DY0 - 1/R0*DY0
  }
  
  if(!doing.CV & verbose) print(paste0('Solve EIF: ', mean(DY) ))

  
  # estimated variance for txt specific means or if break the match	
  var.R1 <- var(DY1) /J
  var.R0 <- var(DY0) / J
  var.break <- var(DY) /J
  
  if( 'pair' %in% colnames(Vdata) ){
    # estimated variance if preserve the match
    pairC <- aggregate(Vdata, by=list(Vdata$id), mean)[,'pair']
    pairs <- unique(pairC)
    n.pairs <- length(pairs)
    DY.paired <-  rep(NA, n.pairs)
    for(i in 1:n.pairs){		
      these<- pairC %in% pairs[i] 
      DY.paired[i]<- 0.5*sum(DY[ these] )			
    }
    
    var.pair <- var(DY.paired) / n.pairs
  } else{
    DY.paired <- var.pair <- NA
  }

  
  
  list(R1=R1, R0=R0, DY1=DY1, var.R1=var.R1, DY0=DY0, var.R0=var.R0, 
       DY=DY, var.break=var.break, 
       DY.paired=DY.paired, var.pair=var.pair)
}





#-----------------------------------------------------#-----------------------------------------------------
# get.inference: function to calculate two-sided confidence intervals
#     & test the null hypothesis with a one-sided test
#	input: 
#		goal (aRR= arithmetic risk ratio; otherwise RD)
#   psi (true value)
#   psi.hat (estimate)
#   se (standard error)
#		df (degrees of freedom if using a Student's t-dist ) 
#		sig.level (significance level)
#   one.sided (if one-sided test)
#   alt.smaller (if under the alternative hypothesis, the intervention reduces the mean outcome)
# output: 
#   point estimate, lower 95%CI, upper 95%CI, standard error estimate, pvalue,
#   if in a simulation (truth psi is known), bias, coverage
#   reject: indicator that reject the null 
# 	note: if goal={aRR,OR} se is on log-scale
#-----------------------------------------------------#-----------------------------------------------------	

get.inference <- function(goal='RD', psi=NA, psi.hat, se, df=99, sig.level=0.05, 
                          one.sided=F, alt.smaller=NULL){
  
  # if doing a one-sided test, need to specify the alternative
  # alt.smaller=T if intervention reduces mean outcome
  # alt.smaller=F if intervention increases mean outcome
  if(one.sided & is.null(alt.smaller)){
    print('*****ERROR: For one-sided test, need to specify the direction of the hypo')
  }
  
  # test statistic (on the log-transformed scale if goal= aRR or OR )
  tstat <- psi.hat/se
  
  if(df>40){
    # assume normal distribution
    cutoff <- qnorm(sig.level/2, lower.tail=F)
    # one.sided hypothesis test 
    if(one.sided){
      pval<- pnorm(tstat, lower.tail=alt.smaller) 
    } else{
      pval<- 2*pnorm(abs(tstat), lower.tail=F) 
    }
  }else{
    # use Student's t-distribution
    # print('Using t-distribution')
    cutoff <- qt(sig.level/2, df=df, lower.tail=F)
    # one.sided hypothesis test 
    if(one.sided){
      pval <- pt(tstat, df=df, lower.tail= alt.smaller ) 
    } else{
      pval <- 2*pt(abs(tstat), df=df, lower.tail=F)
    }
  }
  
  # 95% confidence interval 
  CI.lo <- (psi.hat - cutoff*se)
  CI.hi <- (psi.hat + cutoff*se)
  
  # transform back 
  if(goal!='RD'){
    psi.hat<- exp(psi.hat)
    CI.lo <- exp(CI.lo)
    CI.hi <- exp(CI.hi)
  }  
  
  # bias
  bias <- (psi.hat - psi)
  
  # confidence interval coverage
  cover<- ( CI.lo <= psi & psi <= CI.hi )
  # reject the null
  reject <- as.numeric( pval < sig.level  )
  
  data.frame(est=psi.hat,  CI.lo, CI.hi, se=se,  pval, bias, cover, reject)
  
}

