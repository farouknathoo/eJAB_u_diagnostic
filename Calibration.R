# Calibration.R — Test the new adaptive C*(alpha) calibration plot
#
# Runs the prototype from calibration_new.R on the 500-test simulation data.
# Produces:
#   1. A 3-page PDF (calibration curve, C*(alpha) plot, QQ-plot)
#   2. An Analysis.md summary with detection results and ground truth comparison
#
# Run from the RPackage/ directory:
#   Rscript Temp/Calibration.R

# Source the package code (for ejab01, diagnostic_U, etc.) and the prototype
source("Package/ejabT1E/R/ejab_t1e.R")
source("Temp/calibration_new.R")

# --- Settings ---
# alpha: the pre-chosen significance level for T1E detection (Plot 3)
# up:    the upper bound u_p; only results with p < up are used
#        Target proportion at any alpha is alpha/up.
#        e.g. alpha=0.05, up=0.1 => target=0.5 (very high)
#             alpha=0.01, up=0.1 => target=0.1 (more reasonable)
alpha <- 0.05
up    <- 0.10

# --- Load and prepare data ---
sim <- read.csv("Simulations/ejab_simulation_data_500.csv")
sim <- sim[sim$n > 1 & sim$p_value > 0 & sim$p_value < 1 & sim$q >= 1, ]
df_input <- data.frame(p = sim$p_value, n = sim$n, q = sim$q)

# Compute eJAB01 Bayes factors for all results
ejab_vals <- ejab01(df_input$p, df_input$n, df_input$q)

# --- Run new calibration and produce plots ---
# grid_range = c(0, 1): search C values from 0 to 1
# grid_n = 200: number of C candidate values
# n_alpha = 200: number of alpha grid points for Plots 1-2
pdf("Temp/calibration_new_output.pdf", width = 8, height = 6)

result <- calibration_plot_new(df_input$p, ejab_vals, up = up, alpha = alpha,
                                grid_range = c(0, 1), grid_n = 200,
                                n_alpha = 200,
                                n = df_input$n, q = df_input$q)

dev.off()

# --- Detect candidate T1Es using C*(alpha) ---
Cstar <- result$Cstar_at_alpha
predicted_idx <- detect_type1(df_input$p, ejab_vals, alpha, Cstar)

# --- Ground truth comparison (simulation has is_null column) ---
sim$predicted_t1e <- FALSE
sim$predicted_t1e[predicted_idx] <- TRUE
real_t1e <- sim$is_null == 1 & sim$p_value < alpha
n_real       <- sum(real_t1e)
n_predicted  <- sum(sim$predicted_t1e)
n_overlap    <- sum(real_t1e & sim$predicted_t1e)
n_total      <- nrow(sim)

# Sensitivity: of real T1Es, how many did we catch?
sensitivity <- ifelse(n_real > 0, n_overlap / n_real, NA)
# Specificity: of non-T1Es, how many did we correctly leave alone?
n_not_t1e <- n_total - n_real
n_true_neg <- n_not_t1e - (n_predicted - n_overlap)
specificity <- ifelse(n_not_t1e > 0, n_true_neg / n_not_t1e, NA)

# --- Write analysis to file ---
out <- file("Temp/Analysis.md", open = "w")

writeLines(c(
  "# eJAB Type I Error Analysis (Adaptive C*(alpha))",
  "",
  "## Settings",
  "",
  sprintf("- Alpha: %g", alpha),
  sprintf("- Up: %g", up),
  sprintf("- Studies: %d", n_total),
  sprintf("- C grid: [0, 1] with 200 points"),
  "",
  "## C*(alpha) Estimation",
  "",
  sprintf("- C*(alpha = %g): %.4f", alpha, Cstar),
  sprintf("- Method: pointwise grid search minimising (proportion - alpha/up)^2"),
  "",
  "## Detection Results",
  "",
  sprintf("- Predicted T1Es: %d / %d (%.1f%%)",
          n_predicted, n_total, 100 * n_predicted / n_total),
  "",
  "## Ground Truth Comparison",
  "",
  sprintf("- True T1Es (null & p < alpha): %d / %d (%.1f%%)",
          n_real, n_total, 100 * n_real / n_total),
  sprintf("- Correctly predicted: %d / %d (%.1f%%)",
          n_overlap, n_real,
          ifelse(n_real > 0, 100 * n_overlap / n_real, 0)),
  sprintf("- False positives: %d", n_predicted - n_overlap),
  sprintf("- Sensitivity: %.3f", sensitivity),
  sprintf("- Specificity: %.3f", specificity),
  "",
  "## Plots",
  "",
  "- `calibration_new_output.pdf` page 1 -- Calibration curve using adaptive C*(alpha)",
  "- `calibration_new_output.pdf` page 2 -- C*(alpha) vs alpha",
  "- `calibration_new_output.pdf` page 3 -- Diagnostic QQ-plot with best-fit line",
  "",
  "## Candidates",
  "",
  sprintf("See `candidates.csv` (%d candidates)", n_predicted)
), out)

close(out)

# --- Save candidates to CSV ---
if (n_predicted > 0) {
  candidates <- data.frame(
    Index = predicted_idx,
    p_value = df_input$p[predicted_idx],
    eJAB01 = ejab_vals[predicted_idx]
  )
  write.csv(candidates, "Temp/candidates.csv", row.names = FALSE)
}

# --- Console summary ---
cat("C*(alpha =", alpha, ") =", round(Cstar, 4), "\n")
cat("Predicted T1Es:", n_predicted, "\n")
cat("Sensitivity:", round(sensitivity, 3), "\n")
cat("Specificity:", round(specificity, 3), "\n")
cat("Output saved to Temp/calibration_new_output.pdf and Temp/Analysis.md\n")
