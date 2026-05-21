# Bucher indirect treatment comparison for grade >=3 treatment-related adverse events
# ICI + BCG regimens in BCG-naive high-risk NMIBC

setwd("~/NMA Overall survival Kidney Immuno/BladderNMA/Github R scripts")

# ------------------------------------------------------------
# Safety input data
# Event counts are grade >=3 treatment-related adverse events
# from each trial safety population.
# ------------------------------------------------------------

safety_input <- data.frame(
  study = c("ALBAN", "CREST", "POTOMAC"),
  treatment = c("Atezolizumab+BCG", "Sasanlimab+BCG", "Durvalumab+BCG"),
  comparator = c("BCG", "BCG", "BCG"),
  events_treatment = c(58, 102, 71),
  n_treatment = c(255, 350, 336),
  events_comparator = c(22, 22, 13),
  n_comparator = c(250, 349, 339)
)

write.csv(safety_input, "safety_TRAE_input_table.csv", row.names = FALSE)

# ------------------------------------------------------------
# Calculate odds ratios versus BCG for each active regimen
# ------------------------------------------------------------

safety_input$non_events_treatment <- safety_input$n_treatment - safety_input$events_treatment
safety_input$non_events_comparator <- safety_input$n_comparator - safety_input$events_comparator

safety_input$OR_vs_BCG <- with(
  safety_input,
  (events_treatment / non_events_treatment) /
    (events_comparator / non_events_comparator)
)

safety_input$logOR <- log(safety_input$OR_vs_BCG)

safety_input$SE_logOR <- with(
  safety_input,
  sqrt(
    1 / events_treatment +
      1 / non_events_treatment +
      1 / events_comparator +
      1 / non_events_comparator
  )
)

safety_effects <- safety_input[, c(
  "study",
  "treatment",
  "comparator",
  "events_treatment",
  "n_treatment",
  "events_comparator",
  "n_comparator",
  "OR_vs_BCG",
  "logOR",
  "SE_logOR"
)]

write.csv(safety_effects, "safety_OR_vs_BCG_results.csv", row.names = FALSE)

# ------------------------------------------------------------
# Bucher indirect comparison function
# Compares treatment A vs treatment B through common comparator BCG
# ------------------------------------------------------------

bucher_compare <- function(logOR_A, SE_A, logOR_B, SE_B, label_A, label_B) {
  logOR_indirect <- logOR_A - logOR_B
  SE_indirect <- sqrt(SE_A^2 + SE_B^2)
  
  OR <- exp(logOR_indirect)
  lower_CI <- exp(logOR_indirect - 1.96 * SE_indirect)
  upper_CI <- exp(logOR_indirect + 1.96 * SE_indirect)
  z <- logOR_indirect / SE_indirect
  p_value <- 2 * (1 - pnorm(abs(z)))
  
  data.frame(
    comparison = paste(label_A, "vs", label_B),
    OR = OR,
    lower_CI = lower_CI,
    upper_CI = upper_CI,
    p_value = p_value
  )
}

# Extract active regimen effects versus BCG

atezo <- subset(safety_input, treatment == "Atezolizumab+BCG")
sasa <- subset(safety_input, treatment == "Sasanlimab+BCG")
durva <- subset(safety_input, treatment == "Durvalumab+BCG")

# ------------------------------------------------------------
# Pairwise Bucher indirect comparisons between active regimens
# ------------------------------------------------------------

bucher_results <- rbind(
  bucher_compare(
    atezo$logOR, atezo$SE_logOR,
    sasa$logOR, sasa$SE_logOR,
    "Atezolizumab+BCG", "Sasanlimab+BCG"
  ),
  bucher_compare(
    atezo$logOR, atezo$SE_logOR,
    durva$logOR, durva$SE_logOR,
    "Atezolizumab+BCG", "Durvalumab+BCG"
  ),
  bucher_compare(
    sasa$logOR, sasa$SE_logOR,
    durva$logOR, durva$SE_logOR,
    "Sasanlimab+BCG", "Durvalumab+BCG"
  )
)

write.csv(bucher_results, "bucher_safety_itc_results.csv", row.names = FALSE)

# Print outputs

safety_input
safety_effects
bucher_results
