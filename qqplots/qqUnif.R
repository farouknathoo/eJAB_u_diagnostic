#!/usr/bin/env Rscript
#
# Generate QQ-uniform plot
#

library(readr)
library(dplyr)
library(gap)

# ============================================================================
# Configuration
# ============================================================================
ALPHA <- 0.05
DATA_FILE <- "data/ejab_simulation_data.csv"
OUTPUT_DIR <- "figs"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ============================================================================
# Load data and compute u-values (no selection bias)
# ============================================================================
df <- read_csv(DATA_FILE, show_col_types = FALSE)

flags <- df %>%
  filter(p_value <= ALPHA & ejab_value > 1) %>%
  mutate(
    L = 1 - pchisq((n^(1/q) / (n^(1/q) - 1)) * log(n), df = q),
    u = pmin(pmax((p_value - L) / (ALPHA - L), .Machine$double.eps), 
             1 - .Machine$double.eps)
  ) %>%
  arrange(u)

cat("Total flags:", nrow(flags), 
    "| H₀:", sum(flags$is_null == 1),
    "| H₁:", sum(flags$is_null == 0), "\n")

# ============================================================================
# Generate QQ plot
# ============================================================================
png(file.path(OUTPUT_DIR, "qq_uniform.png"), width = 1000, height = 800, res = 120)

qqunif(flags$u,
       col = ifelse(flags$is_null == 1, "steelblue", "firebrick"),
       pch = 20,
       cex = 0.6,
       main = sprintf("QQ-uniform: All Flagged Points (α = %.2f)", ALPHA))

legend("topleft",
       legend = c("True Type I (H₀)", "False flags (H₁)"),
       col = c("steelblue", "firebrick"),
       pch = 20,
       bty = "n")

dev.off()
cat("Saved to:", file.path(OUTPUT_DIR, "qq_uniform.png"), "\n")