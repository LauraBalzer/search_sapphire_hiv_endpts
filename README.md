# search_sapphire_hiv_endpts
This repository provides the code to replicate the HIV-related endpoints for the SEARCH SAPPHIRE trial designed to evaluate the effectivness of a Community Precision Health intervention in rural Kenya and Uganda (NCT05768763).


The primary endpoint was HIV incidence and measured using a recent infection testing algorithm (RITA). The relevant code is in `incidence_run.R` and `incidence_calcs.R`


The secondary HIV-related endpoints included in this repo are
- HIV care cascade & population-level viral suppression amomg all persons with HIV (PHW)
- Biomedical HIV prevention coverage
- Proporton of all adults (regardless of HIV status) with HIV RNA >1000 c/mL & HIV RNA >99,999 c/mL
  
The relevant code is in `secondary_hiv_run.R` and `secondary_hiv_functions.R`. 


This repo also includes code to 
- run simulations evaluate to RITA performance via `simulate_primary_endpoint.R` with results in `simulate_primary_endpoint.csv` (SAP version) and via `simulation_tables.R` with results in `table_I.csv` and `table_II.csv` (updated - Nature version)
- use treatment-blind simulations to confirm Type-I error control for the primary endpoint via `txt_blind_sims.R` with results in `inc_txt_blind_reps5000.rdata`


Helper functions to evaluate effectiveness: `stage2.R`, `tmle.R`, & `aps.R`

The full Statistical Analysis Plan (SAP) is available in `SEARCH Community Precision Health SAP_v2.1.pdf`.
