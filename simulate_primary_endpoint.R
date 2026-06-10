# --- Setup ---
use_seed <- T
base_seed <- 12349

library(data.table)
library(inctools)
library(parallel)
source('stage2.R')
source('tmle.R')

# Design & RITA parameters
K_per_arm <- 8
BigT      <- 730
pH        <- 0.10
n_tested_mu <- 5000
niter     <- 1000

# These are required for incprops but don't affect the point estimate 
RSE_MDRI <- 0.139
RSE_FRR <- 0.961
RSE_PrevH <- 0.1
RSE_PrevR <- 0.3


# Helpers
sim_comm <- function(I_true_arm, n_tested_mu, pH, MDRI, FRR) {
  I_true <- I_true_arm * rnorm(1, mean = 1, sd = I_comm_sd)  
  n_tested <- max(rpois(1, n_tested_mu), 1L)
  n_hivpos <- rbinom(1, n_tested, pH)
  PrevH_obs <- n_hivpos / n_tested
  # Consistent with incprops (Kassanjee):
  base <- if (PrevH_obs > 0) (1 - PrevH_obs) / PrevH_obs else Inf
  adj <- MDRI/365 - FRR * BigT/365
  p_recent <- FRR + I_true * adj * base
  p_recent <- min(max(p_recent, 0), 0.999)
  n_recent <- if (n_hivpos > 0) rbinom(1, n_hivpos, p_recent) else 0L
  data.table(n_tested, n_hivpos, n_recent)
}

inc_point <- function(PrevH, PrevR, MDRI, FRR, BigT) {
  if (PrevR == 0) {
    #avoid errors in incprops
    return(0)
  }
  out <- incprops(
    PrevH = PrevH, RSE_PrevH = RSE_PrevH,
    PrevR = PrevR, RSE_PrevR = RSE_PrevR,
    MDRI  = MDRI,  RSE_MDRI  = RSE_MDRI,
    FRR   = FRR,   RSE_FRR   = RSE_FRR,
    BigT  = BigT, BS_Count = 0, Boot = FALSE,
    BMest = "same.test"
  )
  inc <- as.numeric(out$Incidence.Statistics$Incidence[1])
  inc <- pmax(inc, 0)
  return(inc)
}

sim_once <- function(iter) {
  if (use_seed) {
    set.seed(1000 * base_seed + iter)
  }
  # 8 C + 8 I communities
  ctl <- rbindlist(lapply(1:K_per_arm, function(i) sim_comm(I_true_C, n_tested_mu, pH, MDRI_true, FRR_true_C)))
  ctl[, arm := "C"]
  ctl[, cluster := .I]
  int <- rbindlist(lapply(1:K_per_arm, function(i) sim_comm(I_true_I, n_tested_mu, pH, MDRI_true, FRR_true_I)))
  int[, arm := "I"]
  int[, cluster := .I]
  dt <- rbind(ctl, int)
  
  # Proportions
  dt[, PrevH := n_hivpos / pmax(n_tested, 1)]
  dt[, PrevR := fifelse(n_hivpos > 0, n_recent / n_hivpos, NA_real_)]
  dt <- dt[is.finite(PrevR)]   # drop clusters with 0 HIV+
  
  dt[arm == 'C', FRR_used := FRR_used_C]
  dt[arm == 'I', FRR_used := FRR_used_I]
  
  # Per-community incidence
  dt[, inc_hat := inc_point(PrevH, PrevR, MDRI_used, FRR_used, BigT), by = .(arm, cluster)]
  df <- as.data.frame(dt[, .(id = 1:16, U = 1, alpha = 1, A = ifelse(arm=='I', 1, 0), Y = inc_hat)])
  
  # for the relative effect -  using unadjusted estimator
  rel <- Stage2(goal='aRR', data.input=df, do.unadjusted=T, 
                psi=I_true_I/I_true_C, one.sided = T, alt.smaller = T)
  
  ratio_est <- rel$est
  ratio_lo <- rel$CI.lo
  ratio_hi <- rel$CI.hi
  
  data.table(
    ratio_est = ratio_est, ratio_lo = ratio_lo, ratio_hi = ratio_hi, ratio_se = rel$se, ratio_reject = rel$reject
  )
  
}

sim <- function() {
  # Run iterations
  res <- rbindlist(mclapply(1:niter, function(i) sim_once(i), mc.cores = 8))

  # True effects (on the target scale)
  true_diff  <- I_true_I - I_true_C
  true_ratio <- I_true_I / I_true_C
  
  # Summarize bias & coverage
  summary_dt <- data.table(
    ratio_bias  = mean(res$ratio_est - true_ratio),
    ratio_sd = sd(res$ratio_est),
    ratio_cov95 = mean(res$ratio_lo <= true_ratio & true_ratio <= res$ratio_hi)
  )
  
  params <- data.table(FRR_true_base,
                       FRR_true_C, 
                       FRR_used_C,
                       FRR_true_I, 
                       FRR_used_I,
                       MDRI_true,
                       MDRI_used,
                       eligible_I,
                       I_comm_sd,
                       I_true_I)
  print(params)
  return(data.table(params, summary_dt))
}

results <- NULL

# Fixed parameters
eligible_C <- 0.07
FRR_used_C <- FRR_used_I <- 0
MDRI_used <- 130
I_true_C  <- 0.005

# Varying parameters

# Note: SAP Table A.1 shows eligible_I = 0.035, I_comm_sd = 0.003
for (MDRI_true in c(130, 198.90)) {
  for (FRR_true_base in c(0, 0.0155)) { #before adjusting for eligible
    for (eligible_I in c(0.035, 0.05, 0.07)) {
      for (I_comm_sd in c(0.001, 0.003)) {
        for (I_true_I  in c(0.0025, 0.00375, 0.005)) {
          FRR_true_C <- FRR_true_base * eligible_C
          FRR_true_I <- FRR_true_base * eligible_I
          results <- rbind(results, sim())
        }
      }
    }
  }
}

print(results)
fwrite(results, file = "simulate_primary_endpoint.csv")
