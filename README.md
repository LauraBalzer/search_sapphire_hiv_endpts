# search_sapphire_hiv_endpts

This repository provides the code to replicate the HIV-related endpoints for the SEARCH SAPPHIRE trial, designed to evaluate the effectiveness of a Community Precision Health intervention in rural Kenya and Uganda (NCT05768763).

The primary endpoint was HIV incidence, measured using a recent infection testing algorithm (RITA). The relevant code is in `incidence_run.R` and `incidence_calcs.R`.

The secondary HIV-related endpoints included in this repository are:

- HIV care cascade and population-level viral suppression among all persons with HIV (PWH)
- Biomedical HIV prevention coverage
- Proportion of all adults (regardless of HIV status) with HIV RNA >1,000 copies/mL and >99,999 copies/mL

The relevant code is in `secondary_hiv_run.R` and `secondary_hiv_functions.R`.

This repository also includes code to:

- evaluate the performance of the primary analysis by simulation:
  - `simulate_primary_endpoint.R`, with results in `simulate_primary_endpoint.csv`: robustness to RITA parameter specifications, as in the Statistical Analysis Plan (SAP);
  - `simulation_tables.R`, with results in `table_I.csv` and `table_II.csv`: updated simulations with pair-matched (restricted) randomization and realistic between-community heterogeneity, analysed with and without breaking the matched pairs (Supplementary Tables S19 and S20);
- run the pre-unblinding treatment-blind simulations, which randomly reassigned arms across communities without regard to the matched pairs: `txt_blind_sims.R`, with results in `inc_txt_blind_reps5000.rdata`.

Helper functions for estimating effects: `stage2.R`, `tmle.R` and `aps.R`.

The full Statistical Analysis Plan is in `SEARCH Community Precision Health SAP_v2.1.pdf`.
