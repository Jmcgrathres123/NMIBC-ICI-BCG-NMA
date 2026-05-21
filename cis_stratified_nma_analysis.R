# CIS-stratified network meta-analysis inputs and P-score ranking
# ICI + BCG regimens in BCG-naive high-risk NMIBC

setwd("~/NMA Overall survival Kidney Immuno/BladderNMA/Github R scripts")

library(netmeta)

# CIS-stratified HR/CI inputs
# HRs and 95% CIs correspond to the trial-reported CIS subgroup analyses used in Figure 4.

cis_input <- data.frame(
  analysis = c(
    "CIS_absent", "CIS_absent", "CIS_absent",
    "CIS_present", "CIS_present", "CIS_present"
  ),
  study = c(
    "ALBAN", "CREST", "POTOMAC",
    "ALBAN", "CREST", "POTOMAC"
  ),
  treat1 = c(
    "Atezolizumab+BCG", "Sasanlimab+BCG", "Durvalumab+BCG",
    "Atezolizumab+BCG", "Sasanlimab+BCG", "Durvalumab+BCG"
  ),
  treat2 = c("BCG", "BCG", "BCG", "BCG", "BCG", "BCG"),
  HR = c(
    1.12, 0.76, 0.53,
    0.74, 0.53, 1.01
  ),
  lower_CI = c(
    0.75, 0.52, 0.35,
    0.43, 0.29, 0.62
  ),
  upper_CI = c(
    1.68, 1.11, 0.81,
    1.28, 0.97, 1.64
  )
)

# Export CIS-stratified input table for repository
write.csv(cis_input, "cis_stratified_input_hr_ci_table.csv", row.names = FALSE)

# Convert HRs and CIs to log HR and standard error
cis_input$TE <- log(cis_input$HR)
cis_input$seTE <- (log(cis_input$upper_CI) - log(cis_input$lower_CI)) / (2 * 1.96)

# CIS-absent analysis
cis_absent <- subset(cis_input, analysis == "CIS_absent")

nma_cis_absent <- netmeta(
  TE = TE,
  seTE = seTE,
  treat1 = treat1,
  treat2 = treat2,
  studlab = study,
  data = cis_absent,
  sm = "HR",
  random = TRUE,
  common = FALSE,
  reference.group = "BCG"
)

summary(nma_cis_absent)

rank_cis_absent <- netrank(nma_cis_absent, small.values = "good")
rank_cis_absent

pscore_cis_absent <- data.frame(
  analysis = "CIS_absent",
  Treatment = names(rank_cis_absent$Pscore.random),
  P_score = as.numeric(rank_cis_absent$Pscore.random)
)

# CIS-present analysis
cis_present <- subset(cis_input, analysis == "CIS_present")

nma_cis_present <- netmeta(
  TE = TE,
  seTE = seTE,
  treat1 = treat1,
  treat2 = treat2,
  studlab = study,
  data = cis_present,
  sm = "HR",
  random = TRUE,
  common = FALSE,
  reference.group = "BCG"
)

summary(nma_cis_present)

rank_cis_present <- netrank(nma_cis_present, small.values = "good")
rank_cis_present

pscore_cis_present <- data.frame(
  analysis = "CIS_present",
  Treatment = names(rank_cis_present$Pscore.random),
  P_score = as.numeric(rank_cis_present$Pscore.random)
)

# Export CIS-stratified P-score results
pscore_cis_stratified <- rbind(pscore_cis_absent, pscore_cis_present)

write.csv(pscore_cis_stratified, "cis_stratified_pscore_results.csv", row.names = FALSE)

pscore_cis_stratified

