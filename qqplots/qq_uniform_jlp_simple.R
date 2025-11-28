#!/usr/bin/env Rscript

# QQ-uniform plot for JLP05, restricted to p < 0.05 and JAB > 1.
# Uses the JAB>1-derived lower bound L(N,k) to compute selection-adjusted u-values.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

ALPHA        <- 0.05
DATA_FILE    <- "data/ejabclinicaltrials/csv/JLP05.csv"
OUTPUT_DIR   <- "figs"
NEW_DATA_DIR <- "data/new"

if (!dir.exists(OUTPUT_DIR))   dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(NEW_DATA_DIR)) dir.create(NEW_DATA_DIR, recursive = TRUE)

# -------------------------------------------------------------------
# Helper: infer k (df) from model, I, R, C, matching JAB01
# -------------------------------------------------------------------
infer_k <- function(model, I, R, C) {
  out <- rep(NA_integer_, length(model))

  is_t1 <- model %in% c(
    "t-test", "linear_regression",
    "logistic_regression", "cox",
    "wilcoxon", "mann_whitney"
  )
  out[is_t1] <- 1L

  is_anova_like <- model %in% c("anova", "kruskal_wallis", "repeated_measures")
  out[is_anova_like] <- ifelse(
    !is.na(I[is_anova_like]),
    I[is_anova_like] - 1L,
    NA_integer_
  )

  is_chisq <- model == "chi_squared"
  out[is_chisq] <- ifelse(
    !is.na(R[is_chisq]) & !is.na(C[is_chisq]),
    (R[is_chisq] - 1L) * (C[is_chisq] - 1L),
    NA_integer_
  )

  out
}

# -------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------
df <- read_csv(DATA_FILE, show_col_types = FALSE)

needed_cols <- c("pValue", "JAB", "N", "model", "I", "R", "C")
missing_cols <- setdiff(needed_cols, names(df))
if (length(missing_cols) > 0L) {
  stop("Missing expected columns in JLP05: ",
       paste(missing_cols, collapse = ", "))
}

# -------------------------------------------------------------------
# Restrict to p < 0.05 and JAB > 1, compute u-values
# -------------------------------------------------------------------
df2 <- df %>%
  mutate(k = infer_k(model, I, R, C))

if (any(is.na(df2$k))) {
  warning("Some rows have NA k (cannot infer df); they will be dropped.")
}

flags <- df2 %>%
  filter(
    !is.na(pValue),
    !is.na(JAB),
    !is.na(N),
    !is.na(k),
    pValue < ALPHA,
    JAB > 1
  ) %>%
  mutate(
    A = (N^(1 / k) - 1) / (N^(1 / k)),
    Q_thresh = log(N) / A,
    L = 1 - pchisq(Q_thresh, df = k),

    u_raw = (pValue - L) / (ALPHA - L),
    u = pmin(pmax(u_raw, .Machine$double.eps), 1 - .Machine$double.eps)
  ) %>%
  arrange(u)

cat("Total points with p < ", ALPHA, " and JAB > 1: ",
    nrow(flags), "\n", sep = "")

if (nrow(flags) == 0L) {
  stop("No rows satisfy p < 0.05 and JAB > 1.")
}

# save the subset
out_file <- file.path(NEW_DATA_DIR, "JLP05_p_lt_0.05_JAB_gt_1_u.csv")
write_csv(flags, out_file)
cat("Flagged data with u-values written to: ", out_file, "\n", sep = "")

# -------------------------------------------------------------------
# QQ-uniform plot
# -------------------------------------------------------------------
m <- nrow(flags)
expected <- (seq_len(m) - 0.5) / m
observed <- flags$u

png(file.path(OUTPUT_DIR, "qq_uniform_JLP05_p_lt_0.05_JAB_gt_1.png"),
    width = 1000, height = 800, res = 120)

plot(expected, observed,
     xlab = "Expected U(0,1) quantiles",
     ylab = "Observed u-values",
     main = "QQ-uniform (JLP05): p < 0.05 & JAB > 1",
     pch  = 20,
     cex  = 0.6)

abline(0, 1, lty = 2)

dev.off()

cat("Saved QQ plot to: ",
    file.path(OUTPUT_DIR, "qq_uniform_JLP05_p_lt_0.05_JAB_gt_1.png"),
    "\n", sep = "")
