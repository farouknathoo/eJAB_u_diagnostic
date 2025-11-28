#!/usr/bin/env Rscript
#
# QQ-uniform plot for already-flagged eJAB data (simulation or CTG),
# assuming the CSV already contains 'u' (and optionally 'B_u').
#

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ============================================================================
# Parse command-line arguments
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)

type <- NULL          # "simulation" or "ctg"
data_file <- NULL
output_dir <- "figs"
alpha_label <- "0.05"  # just for the title

for (a in args) {
  if (a == "--simulation") {
    type <- "simulation"
  } else if (a == "--ctg") {
    type <- "ctg"
  } else if (grepl("^--data-file=", a)) {
    data_file <- sub("^--data-file=", "", a)
  } else if (grepl("^--output-dir=", a)) {
    output_dir <- sub("^--output-dir=", "", a)
  } else if (grepl("^--alpha=", a)) {
    alpha_label <- sub("^--alpha=", "", a)
  }
}

if (is.null(type)) {
  stop("You must specify one of --simulation or --ctg", call. = FALSE)
}

if (is.null(data_file)) {
  data_file <- if (type == "simulation") {
    "data/ejab_simulation_flagged_with_u_B.csv"
  } else {
    "data/CTG_clean_flagged_with_u_B.csv"
  }
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================================
# Load data (already flagged) and prepare u-values
# ============================================================================
df <- read_csv(data_file, show_col_types = FALSE)

if (!"u" %in% names(df)) {
  stop("Input file must contain a 'u' column (precomputed).", call. = FALSE)
}

# Keep non-missing u and sort
df <- df %>%
  filter(!is.na(u)) %>%
  arrange(u)

n <- nrow(df)

if ("is_null" %in% names(df)) {
  cat("Rows:", n,
      "| H0:", sum(df$is_null == 1, na.rm = TRUE),
      "| H1:", sum(df$is_null == 0, na.rm = TRUE), "\n")
} else {
  cat("Rows:", n, "\n")
}

if (n == 0) {
  warning("No rows with non-missing u; no plot produced.")
  quit(status = 0)
}

# Clip u to (0,1) just in case
df <- df %>%
  mutate(
    u = pmin(pmax(u, .Machine$double.eps),
             1 - .Machine$double.eps)
  )

# ============================================================================
# Compute theoretical quantiles and make QQ plot on actual [0,1] scale
# ============================================================================
theo <- (seq_len(n) - 0.5) / n
obs  <- df$u

out_file <- if (type == "simulation") {
  file.path(output_dir, "qq_uniform_simulation.png")
} else {
  file.path(output_dir, "qq_uniform_ctg.png")
}

png(out_file, width = 1000, height = 800, res = 120)

if ("is_null" %in% names(df) && any(!is.na(df$is_null))) {
  point_col <- ifelse(df$is_null == 1, "steelblue", "firebrick")
} else {
  point_col <- "steelblue"
}

plot(
  theo, obs,
  xlab = "Theoretical quantiles (Uniform(0,1))",
  ylab = "Observed u",
  xlim = c(0, 1),
  ylim = c(0, 1),
  pch  = 20,
  col  = point_col,
  main = sprintf(
    "QQ-uniform (%s, α = %s)",
    ifelse(type == "simulation", "Simulation", "CTG"),
    alpha_label
  )
)

abline(0, 1, lty = 2, col = "gray40")

if ("is_null" %in% names(df) && any(!is.na(df$is_null))) {
  legend("topleft",
         legend = c("True Type I (H0)", "False flags (H1)"),
         col    = c("steelblue", "firebrick"),
         pch    = 20,
         bty    = "n")
}

dev.off()
cat("Saved to:", out_file, "\n")
