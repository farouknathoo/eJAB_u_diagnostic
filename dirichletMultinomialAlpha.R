library(gtools)
source("CTGJAB.R")
source("rfuncs/DirichletMultinomial.R")

# Define the range of alpha values
alpha_values <- c(
  0.0001, 0.001, 0.005,
  0.01, 0.025, 0.05, 0.1, 0.25,
  0.5, 0.75, 0.9
)

# Create a directory for saving results if it doesn't exist
if (!dir.exists("csv")) {
  dir.create("csv")
}

# Initialize a data frame to store all results
results_df <- data.frame(
  alpha = numeric(),
  posterior_mean = numeric(),
  posterior_se = numeric(),
  stringsAsFactors = FALSE
)

# Run sequentially for each alpha value
for (alpha in alpha_values) {
  # Run posterior analysis for the current alpha value
  posterior_results <- posterior_analysis(study3_summary, 5000, alpha)

  # Add the results to the data frame
  results_df <- rbind(results_df, data.frame(
    alpha = alpha,
    posterior_mean = posterior_results$mean,
    posterior_se = posterior_results$se
  ))

  # Print progress
  cat("Completed analysis for alpha =", alpha, "\n")
}

# Write all results to a single CSV file
write.csv(results_df, "csv/posterior_results.csv", row.names = FALSE)
