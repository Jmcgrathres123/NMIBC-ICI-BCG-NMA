# Network meta-analysis of ICI + BCG regimens in BCG-naive high-risk NMIBC
# Primary NMA and P-score ranking reproducibility script

setwd("~/NMA Overall survival Kidney Immuno/BladderNMA/Github R scripts")

library(netmeta)
# Input data
# HRs and 95% CIs correspond to the values used in Figure 2 of the manuscript.

nma_input <- data.frame(
  study = c("ALBAN", "CREST", "POTOMAC"),
  treat1 = c("Atezolizumab+BCG", "Sasanlimab+BCG", "Durvalumab+BCG"),
  treat2 = c("BCG", "BCG", "BCG"),
  HR = c(0.98, 0.68, 0.68),
  lower_CI = c(0.71, 0.49, 0.50),
  upper_CI = c(1.36, 0.94, 0.93)
)

# Export input table for repository
write.csv(nma_input, "nma_input_hr_ci_table.csv", row.names = FALSE)
# Convert HRs and 95% CIs to log HR and standard error

nma_input$TE <- log(nma_input$HR)
nma_input$seTE <- (log(nma_input$upper_CI) - log(nma_input$lower_CI)) / (2 * 1.96)
# Random-effects network meta-analysis
# BCG is the reference comparator.

nma.re <- netmeta(
  TE = TE,
  seTE = seTE,
  treat1 = treat1,
  treat2 = treat2,
  studlab = study,
  data = nma_input,
  sm = "HR",
  random = TRUE,
  fixed = FALSE,
  reference.group = "BCG"
)

summary(nma.re)
# P-score treatment ranking
# Lower HR indicates greater benefit, so small.values = "good".

rank.re <- netrank(nma.re, small.values = "good")
rank.re

# Export P-score table

pscore_table <- data.frame(
  Treatment = names(rank.re$Pscore.random),
  P_score = as.numeric(rank.re$Pscore.random)
)

write.csv(pscore_table, "pscore_results.csv", row.names = FALSE)

pscore_table
