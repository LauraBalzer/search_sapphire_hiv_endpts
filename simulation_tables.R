# =============================================================================
# simulation_tables.R
#
# SEARCH Community Precision Health trial (NCT05768763): simulation study of
# the primary analysis of HIV incidence (Two-Stage approach with recency-based
# community incidence estimates). Reproduces the two simulation tables in the
# Supplementary Information:
#
#   Table I   Restricted (pair-matched) randomization: coverage, type I error,
#             power and CI width of the relative effect, analysing each
#             simulated trial with matched pairs broken (pre-specified primary
#             analysis) and maintained, at three levels of matching quality.
#   Table II  Robustness of the primary analysis to misspecification of the
#             recent infection testing algorithm (RITA) parameters (mean
#             duration of recent infection [MDRI], false recency rate [FRR]);
#             extends Table A.1 of the statistical analysis plan.
#
# Data-generating process (one simulated trial)
#   - 8 matched pairs of communities (16 communities); within each pair one
#     community is randomized to intervention.
#   - HIV prevalence: pair-level value ~ Uniform(4%, 17%); each community varies
#     around its pair's value (lognormal, coefficient of variation 0.15).
#   - True HIV incidence = arm incidence x pair effect x community effect, both
#     lognormal with mean 1. Their combined coefficient of variation is 0.4;
#     rho is the share of the log-scale variance shared within pairs (matching
#     quality: 0.2 low, 0.5 moderate, 0.8 high). A CV of 0.4 reproduces the
#     variation in estimated incidence observed across control communities once
#     sampling variability is accounted for.
#   - ~5,000 adults tested per community (Poisson); persons with HIV ~
#     Binomial(tested, prevalence); recent infections ~ Binomial(persons with
#     HIV, Pr[recent]), with Pr[recent] = FRR + incidence x (MDRI - FRR x T) x
#     (1 - prevalence) / prevalence, T = 730 days. The true FRR applies to
#     persons eligible for recency testing (FRR x eligible fraction).
#
# Analysis of each simulated trial
#   - Stage 1: community incidence by the Kassanjee estimator (inctools::incprops)
#     with MDRI 130 days and FRR 0, as in the primary analysis.
#   - Stage 2: unadjusted comparison of mean community incidence between arms
#     (relative effect), influence-curve standard errors; pairs broken (t, 14
#     df) or maintained (t, 7 df); one-sided test at alpha 0.05.
#
# Requirements: R >= 4.1; packages data.table, inctools; stage2.R and tmle.R
# from the trial's analysis code (set code_dir below). Uses parallel::mclapply
# (forking; on Windows set N_CORES <- 1). Each simulated trial has its own
# seed, so results do not depend on the number of cores. 1,000 trials per
# scenario, 18 scenarios
#
# Output: table_I.csv and table_II.csv in the working directory, also printed.
# =============================================================================

library(data.table)
library(inctools)
library(parallel)

code_dir <- "."                          # folder containing stage2.R and tmle.R
source(file.path(code_dir, "stage2.R"))
source(file.path(code_dir, "tmle.R"))

# ---- settings -----------------------------------------------------------------
N_ITER       <- 1000                     # simulated trials per scenario
N_CORES      <- 8
BASE_SEED    <- 30117
N_PAIRS      <- 8
N_TESTED     <- 5000                     # mean adults tested per community
BIG_T        <- 730                      # RITA post-infection time cutoff (days)
INC_CV       <- 0.4                      # between-community CV of true incidence
PH_WITHIN_CV <- 0.15                     # within-pair CV of HIV prevalence
MDRI_USED    <- 130                      # RITA parameters used in the analysis
FRR_USED     <- 0
RSE <- list(MDRI = 0.139, FRR = 0.961, PrevH = 0.1, PrevR = 0.3)   # not used for point estimates

# ---- helpers ------------------------------------------------------------------
# lognormal multiplier with mean 1 and log-scale SD sdlog
LN <- function(n, sdlog) if (sdlog == 0) rep(1, n) else rlnorm(n, -sdlog^2 / 2, sdlog)

IncPoint <- function(PrevH, PrevR) {
  if (PrevR == 0) return(0)
  out <- incprops(PrevH = PrevH, RSE_PrevH = RSE$PrevH, PrevR = PrevR, RSE_PrevR = RSE$PrevR,
                  MDRI = MDRI_USED, RSE_MDRI = RSE$MDRI, FRR = FRR_USED, RSE_FRR = RSE$FRR,
                  BigT = BIG_T, BS_Count = 0, Boot = FALSE, BMest = "same.test")
  pmax(as.numeric(out$Incidence.Statistics$Incidence[1]), 0)
}

# one simulated trial under scenario s (a one-row list of parameters)
SimTrial <- function(s) {
  v <- log(1 + INC_CV^2)
  sd_between <- sqrt(s$rho * v); sd_within <- sqrt((1 - s$rho) * v)
  d <- data.table(pair = 1:N_PAIRS, pH_pair = runif(N_PAIRS, 0.04, 0.17),
                  m_pair = LN(N_PAIRS, sd_between))[rep(1:N_PAIRS, each = 2)]
  d[, A := as.vector(sapply(1:N_PAIRS, function(i) sample(0:1)))]      # within-pair randomization
  d[, m_comm := LN(.N, sd_within)]
  d[, pH := pmin(pH_pair * LN(.N, sqrt(log(1 + PH_WITHIN_CV^2))), 0.5)]
  d[, I_true := fifelse(A == 1, s$inc_I, s$inc_C) * m_pair * m_comm]
  d[, FRR_true := s$frr_true * fifelse(A == 1, s$elig_I, s$elig_C)]
  d[, n_tested := pmax(rpois(.N, N_TESTED), 1L)]
  d[, n_hivpos := rbinom(.N, n_tested, pH)]
  d[, PrevH := n_hivpos / n_tested]
  d[, p_recent := pmin(pmax(FRR_true + I_true * (s$mdri_true / 365 - FRR_true * BIG_T / 365) *
                              (1 - PrevH) / PrevH, 0), 0.999)]
  d[, n_recent := rbinom(.N, n_hivpos, p_recent)]
  d[, inc_hat := IncPoint(PrevH, n_recent / n_hivpos), by = seq_len(nrow(d))]
  d
}

# analyse one trial with pairs broken and maintained
SimOnce <- function(iter, s) {
  set.seed(1000 * BASE_SEED + iter)
  d <- SimTrial(s)
  df <- as.data.frame(d[, .(id = seq_len(.N), pair, U = 1, alpha = 1, A, Y = inc_hat)])
  psi <- s$inc_I / s$inc_C
  fit <- lapply(c(broken = TRUE, maintained = FALSE), function(b)
    Stage2(goal = "aRR", data.input = df, do.unadjusted = TRUE, break.match = b,
           psi = psi, one.sided = TRUE, alt.smaller = TRUE))
  data.table(analysis = names(fit),
             lo = sapply(fit, `[[`, "CI.lo"), hi = sapply(fit, `[[`, "CI.hi"),
             est = sapply(fit, `[[`, "est"), reject = sapply(fit, `[[`, "reject"))
}

RunScenario <- function(s) {
  r <- rbindlist(mclapply(1:N_ITER, SimOnce, s = s, mc.cores = N_CORES))
  psi <- s$inc_I / s$inc_C
  r[is.finite(est) & is.finite(lo) & is.finite(hi),
    .(n_ok = .N, coverage = mean(lo <= psi & psi <= hi), reject = mean(reject),
      ci_width_log = mean(log(hi) - log(lo))), by = analysis][, names(s) := s][]
}

# ---- scenarios ----------------------------------------------------------------
# Table I: incidence and LAg-eligible fractions as observed in the trial
# (control 0.19, intervention 0.055 per 100 person-years; null: 0.19 in both);
# true RITA parameters equal to those used (MDRI 130 days, FRR 0)
grid_I <- CJ(table = "I", rho = c(0.2, 0.5, 0.8), inc_I = c(0.00055, 0.0019), sorted = FALSE)
grid_I[, `:=`(inc_C = 0.0019, elig_C = 0.058, elig_I = 0.041, mdri_true = 130, frr_true = 0)]

# Table II: SAP Table A.1 settings (control 0.5 per 100 person-years), low
# matching quality; true FRR 0 or 0.0155 and true MDRI 130 or 198.9 days
grid_II <- CJ(table = "II", frr_true = c(0, 0.0155), mdri_true = c(130, 198.9),
              inc_I = c(0.0025, 0.00375, 0.005), sorted = FALSE)
grid_II[, `:=`(rho = 0.2, inc_C = 0.005, elig_C = 0.07, elig_I = 0.035)]

grid <- rbind(grid_I, grid_II, use.names = TRUE)
results <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  cat("scenario", i, "of", nrow(grid), "\n")
  RunScenario(as.list(grid[i]))
}))

# ---- tables -------------------------------------------------------------------
F3 <- function(p, n) sprintf("%.3f (%.3f)", p, sqrt(p * (1 - p) / n))   # estimate (Monte Carlo SE)

t1 <- results[table == "I"]
t1[, setting := fifelse(inc_I == inc_C, "null", "effect")]
table_I <- dcast(t1, rho + analysis ~ setting,
                 value.var = c("coverage", "n_ok", "reject", "ci_width_log"))
table_I <- table_I[order(rho, -xtfrm(analysis == "broken"))][, .(
  `Matching quality (rho)` = rho,
  Analysis = fifelse(analysis == "broken", "Pairs broken (pre-specified)", "Pairs maintained"),
  `Coverage, null` = F3(coverage_null, n_ok_null),
  `Type I error` = sprintf("%.3f", reject_null),
  `Coverage, RR 0.29` = F3(coverage_effect, n_ok_effect),
  Power = sprintf("%.3f", reject_effect),
  `Mean CI width (log), RR 0.29` = sprintf("%.2f", ci_width_log_effect))]

t2 <- results[table == "II" & analysis == "broken"]
table_II <- t2[, .(`True FRR` = frr_true, `FRR used` = FRR_USED, `True MDRI` = mdri_true,
                   `MDRI used` = MDRI_USED, `Intervention incidence` = inc_I,
                   `95% CI coverage` = F3(coverage, n_ok))]

cat("\nTable I. Restricted randomization, with and without breaking the match\n")
print(table_I, row.names = FALSE)
cat("\nTable II. Robustness to misspecification of RITA parameters (pairs broken, rho = 0.2)\n")
print(table_II, row.names = FALSE)
cat(sprintf("One-sided type I error (intervention incidence = control): %.3f to %.3f\n",
            t2[inc_I == inc_C, min(reject)], t2[inc_I == inc_C, max(reject)]))

fwrite(table_I, "table_I.csv")
fwrite(table_II, "table_II.csv")
