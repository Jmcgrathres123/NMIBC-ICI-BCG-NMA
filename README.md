# NMIBC-ICI-BCG-NMA

This repository contains reproducibility files for the network meta-analysis of immune checkpoint inhibitor plus BCG regimens in BCG-naïve, high-risk non-muscle-invasive bladder cancer.

## Contents

- `nma_pscore_analysis.R`: Primary network meta-analysis and P-score ranking script.
- `nma_input_hr_ci_table.csv`: Trial-level HR/CI inputs for the primary NMA.
- `pscore_results.csv`: Primary P-score output.
- `cis_stratified_nma_analysis.R`: CIS-stratified subgroup NMA and P-score script.
- `cis_stratified_input_hr_ci_table.csv`: CIS-stratified HR/CI inputs.
- `cis_stratified_pscore_results.csv`: CIS-stratified P-score output.
- `bucher_safety_itc.R`: Bucher indirect treatment comparison for grade ≥3 treatment-related adverse events.
- `safety_TRAE_input_table.csv`: Safety input data for grade ≥3 treatment-related adverse events.
- `safety_OR_vs_BCG_results.csv`: Odds ratios for active regimens versus BCG.
- `bucher_safety_itc_results.csv`: Bucher indirect comparison results.
- `rob2_domain_response.xlsx`: RoB 2 domain-level responses.

## Software

Analyses were performed in R using the `netmeta` package.
