#!/usr/bin/env Rscript

# QQ-uniform plot for JLP05, restricted to p < 0.05 and JAB > 1.
# Uses the JAB>1-derived lower bound L(N,k) to compute selection-adjusted u-values.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

ALPHA        <- 0.05
DATA_FILE    <- "data/JLP05.csv"
OUTPUT_DIR   <- "figs"
NEW_DATA_DIR <- "data"

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
# QQ-uniform plots with trimming by distance from y = x
# -------------------------------------------------------------------
m <- nrow(flags)

# base QQ for full data
expected_full <- (seq_len(m) - 0.5) / m
observed_full <- flags$u
dist_full <- abs(observed_full - expected_full)

flags_plot <- flags %>%
  mutate(
    expected_full = expected_full,
    dist = dist_full
  )

# Trim percentages (of FURTHEST points to drop)
trim_pcts <- c(0, 5, 10, 20, 30, 40, 50)

base_plot_file <- file.path(OUTPUT_DIR, "qq_uniform_JLP05_p_lt_0.05_JAB_gt_1.png")

for (tp in trim_pcts) {
  if (tp <= 0) {
    df_trim <- flags_plot
  } else {
    n_keep <- max(1L, floor(m * (1 - tp / 100)))
    keep_idx <- order(flags_plot$dist)   # sorted by distance (closest first)
    keep_idx <- keep_idx[seq_len(n_keep)]
    df_trim <- flags_plot[keep_idx, , drop = FALSE]
  }

  # sort by u and recompute expected for the trimmed sample
  df_trim <- df_trim %>%
    arrange(u) %>%
    mutate(
      expected = (seq_len(n()) - 0.5) / n(),
      observed = u
    )

  if (tp == 0) {
    out_path  <- base_plot_file
    main_title <- "QQ-uniform (JLP05): p < 0.05 & JAB > 1 (no trimming)"
  } else {
    out_name  <- sprintf("qq_uniform_JLP05_trim_%02d_pct.png", tp)
    out_path  <- file.path(OUTPUT_DIR, out_name)
    main_title <- sprintf("QQ-uniform (JLP05): trimmed %d%% farthest points", tp)
  }

  png(out_path, width = 1000, height = 800, res = 120)
  plot(df_trim$expected, df_trim$observed,
       xlab = "Expected U(0,1) quantiles",
       ylab = "Observed u-values",
       main = main_title,
       pch  = 20,
       cex  = 0.6)
  abline(0, 1, lty = 2)
  dev.off()

  cat("Saved QQ plot (trim ", tp, "%) to: ",
      out_path, "\n", sep = "")
}
