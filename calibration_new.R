# calibration_new.R — Adaptive C*(alpha) calibration plot prototype
#
# New approach to calibration: instead of enforcing a single constant C*
# across all alpha values (which overshoots/undershoots at different alpha),
# we estimate a separate C*(alpha) for each alpha on a grid in [0, u_p].
#
# At each alpha, C*(alpha) is the value of C that makes the observed
# proportion of NHST/Bayes contradictions closest to the theoretical
# Type I error rate alpha/u_p.
#
# Contradiction := (p <= alpha) AND (eJAB01 > C)
# Proportion    := count(contradictions) / count(p < u_p)
# Target        := alpha / u_p
#
# Three outputs:
#   Plot 1: Calibration curve using C*(alpha)
#           x = alpha, y = observed proportion at the optimal C*(alpha)
#           Reference line has slope 1/u_p; curve should track it closely.
#   Plot 2: C*(alpha) vs alpha
#           Shows how the threshold varies with alpha.
#   Plot 3: Diagnostic QQ-plot at a pre-chosen alpha (e.g. 0.05)
#           Uses C*(alpha) from Plot 2 to detect contradictions, computes
#           the diagnostic U_i values, and plots against Unif(0,1).
#           Best-fit line (OLS) replaces the 45-degree line because C* is
#           estimated (not fixed).


# ============================================================================
# estimate_Cstar_alpha()
#
# Core computation: for each alpha on a grid, find C*(alpha) by grid search
# over candidate C values, minimising (observed_proportion - alpha/up)^2.
# ============================================================================

#' @param p         Numeric vector of p-values (all results, not pre-filtered)
#' @param ejab      Numeric vector of eJAB01 values (same length as p)
#' @param up        Upper bound u_p (default 0.1). Only results with p < up
#'                  contribute to the denominator N.
#' @param grid_range c(lower, upper) for the C search grid (default c(0, 1)).
#'                  Must cover the range where meaningful C* values live.
#'                  If too narrow, C*(alpha) may be pinned at a boundary.
#' @param grid_n    Number of C grid points (default 200)
#' @param n_alpha   Number of alpha grid points in (0, up] (default 200)
#' @return list with:
#'   alpha_grid  — the alpha values used
#'   Cstar_alpha — C*(alpha) at each alpha
#'   proportions — the achieved proportion at each alpha (using the optimal C)
estimate_Cstar_alpha <- function(p, ejab, up = 0.1,
                                  grid_range = c(0, 1), grid_n = 200,
                                  n_alpha = 200) {

  # Alpha grid: evenly spaced in (0, up], excluding 0 (target would be 0)
  alpha_grid <- seq(0, up, length.out = n_alpha + 1)[-1]

  # C grid: candidate threshold values to search over
  C_grid <- seq(grid_range[1], grid_range[2], length.out = grid_n)

  # Denominator: number of results with p < up (fixed for all alpha)
  N <- sum(p < up)

  # Output vectors
  Cstar_alpha  <- numeric(length(alpha_grid))
  proportions  <- numeric(length(alpha_grid))

  # For each alpha, find the C that makes the proportion closest to alpha/up
  for (i in seq_along(alpha_grid)) {
    a <- alpha_grid[i]
    target <- a / up   # theoretical T1E rate at this alpha

    # Compute proportion of contradictions at each candidate C
    # proportion(C) = count(p <= alpha AND ejab > C) / N
    # This is a decreasing step function of C
    props <- vapply(C_grid, function(C) {
      sum(p <= a & ejab > C) / N
    }, numeric(1))

    # Pick C minimising squared difference from target
    best <- which.min((props - target)^2)
    Cstar_alpha[i] <- C_grid[best]
    proportions[i] <- props[best]
  }

  list(alpha_grid = alpha_grid,
       Cstar_alpha = Cstar_alpha,
       proportions = proportions)
}


# ============================================================================
# diagnostic_qqplot_fit()
#
# QQ-plot of diagnostic U_i values against Unif(0,1) quantiles.
# Uses OLS best-fit line instead of the 45-degree line, because C* is
# estimated (not fixed at 1), so the diagnostics are a function of C*.
# Includes simultaneous confidence band via MC-calibrated Beta order
# statistics (the MC finds the pointwise level that gives the desired
# simultaneous coverage).
# ============================================================================

#' @param U      Numeric vector of diagnostic U_i values
#' @param alpha  Significance level used for detection (displayed in title)
#' @param Cstar  C* value used for detection (displayed in title)
#' @param band   Logical; add simultaneous confidence band? (default TRUE)
#' @param conf   Confidence level for the simultaneous band (default 0.95)
#' @param B      Number of MC simulations for band calibration (default 10000).
#'               Higher B = tighter calibration but slower.
#' @param seed   Random seed for MC reproducibility (default 1)
#' @param ...    Additional arguments passed to plot()
#' @return Invisibly, list with theoretical, observed, intercept, slope,
#'         lower, upper, outside
diagnostic_qqplot_fit <- function(U, alpha = NULL, Cstar = NULL,
                                   band = TRUE, conf = 0.95,
                                   B = 10000, seed = 1, ...) {
  n <- length(U)
  if (n == 0) {
    message("No candidate T1Es detected; cannot produce QQ-plot.")
    return(invisible(NULL))
  }

  # Theoretical quantiles: expected order statistics from Unif(0,1)
  theoretical <- stats::ppoints(n)

  # Observed quantiles: sorted U values
  observed <- sort(U)

  # --- Best-fit line (OLS) ---
  # Computed first so the confidence band can follow it.
  # Since C* is estimated (not fixed), the U_i may not be exactly Unif(0,1)
  # but rather a shifted/scaled version: U_i ~ a + b * Unif(0,1).
  # The OLS line captures this shift/scale.
  fit_line <- lm(observed ~ theoretical)
  a <- coef(fit_line)[1]  # intercept
  b <- coef(fit_line)[2]  # slope

  # --- Simultaneous confidence band (follows best-fit line) ---
  # Step 1: MC-calibrate the pointwise level p* for Beta(i, n+1-i) bands
  #         such that the simultaneous coverage = conf (e.g. 0.95).
  # Step 2: Construct raw Unif(0,1) bands from Beta quantiles.
  # Step 3: Transform bands through the OLS fit: band_new = a + b * band_raw.
  #         This makes the bands follow the best-fit line rather than the
  #         45-degree line, consistent with using OLS instead of y = x.
  #         The simultaneous coverage is preserved because the linear
  #         transformation is monotone (for b > 0).
  lower <- upper <- outside <- NULL
  if (band && n >= 2) {
    set.seed(seed)
    # Simulate B samples of n Uniform order statistics
    U_sim <- matrix(stats::runif(n * B), nrow = n, ncol = B)
    U_sim <- apply(U_sim, 2, sort)
    i <- seq_len(n)

    # For a given pointwise level p, compute the fraction of simulated
    # samples that fall entirely within the Beta(i, n+1-i) bands
    coverage_hat <- function(p) {
      tail <- (1 - p) / 2
      L <- stats::qbeta(tail, i, n + 1 - i)
      U <- stats::qbeta(1 - tail, i, n + 1 - i)
      # A sample is "inside" if ALL n order statistics are within [L, U]
      mean(colSums(U_sim >= L & U_sim <= U) == n)
    }

    # Find p* such that simultaneous coverage = conf
    f <- function(p) coverage_hat(p) - conf
    if (f(conf) >= 0) {
      # If pointwise level = conf already gives enough simultaneous coverage
      p_star <- conf
    } else {
      # Otherwise, search for the right pointwise level (> conf)
      p_star <- stats::uniroot(f, lower = conf, upper = 0.9999, tol = 1e-4)$root
    }

    # Raw Unif(0,1) bands from Beta order statistics
    tail_star <- (1 - p_star) / 2
    lower_raw <- stats::qbeta(tail_star, i, n + 1 - i)
    upper_raw <- stats::qbeta(1 - tail_star, i, n + 1 - i)

    # Transform bands to follow the best-fit line: y = a + b * x
    lower <- a + b * lower_raw
    upper <- a + b * upper_raw

    # Identify observations falling outside the transformed band
    outside_ord <- which(observed < lower | observed > upper)
    outside <- order(U)[outside_ord]
  }

  # --- Build plot title (include alpha and C* if provided) ---
  main_title <- "Diagnostic QQ-Plot"
  if (!is.null(alpha) && !is.null(Cstar)) {
    main_title <- substitute(
      paste("Diagnostic QQ-Plot at ", alpha == aa, ", ", C^"*" == cs),
      list(aa = alpha, cs = round(Cstar, 4)))
  } else if (!is.null(alpha)) {
    main_title <- bquote("Diagnostic QQ-Plot at " * alpha == .(alpha))
  }

  # --- Plot ---
  plot(theoretical, observed,
       xlab = "Theoretical Unif(0,1) Quantiles",
       ylab = "Observed U_i Quantiles",
       main = main_title,
       pch = 16, cex = 0.6, ...)

  abline(fit_line, col = "red", lwd = 2)

  # Draw transformed confidence band and highlight outliers
  if (band && n >= 2) {
    lines(theoretical, lower, col = "grey50", lty = 2)
    lines(theoretical, upper, col = "grey50", lty = 2)
    if (length(outside_ord) > 0) {
      points(theoretical[outside_ord], observed[outside_ord],
             pch = 16, cex = 0.6, col = "red")
    }
  }

  invisible(list(
    theoretical = theoretical,
    observed = observed,
    intercept = a,
    slope = b,
    lower = lower,
    upper = upper,
    outside = outside
  ))
}


# ============================================================================
# calibration_plot_new()
#
# Main function: runs estimate_Cstar_alpha() and produces all 3 plots.
# Plot 1 and 2 use the full alpha grid; Plot 3 uses the pre-chosen alpha
# (e.g. 0.05) and looks up C*(alpha) from the grid.
# ============================================================================

#' @param p         Numeric vector of p-values
#' @param ejab      Numeric vector of eJAB01 values
#' @param up        Upper bound u_p (default 0.1)
#' @param alpha     Pre-chosen significance level for Plot 3 QQ-plot (default 0.05).
#'                  This is the alpha at which T1Es will be detected.
#' @param grid_range c(lower, upper) for C grid (default c(0, 1))
#' @param grid_n    Number of C grid points (default 200)
#' @param n_alpha   Number of alpha grid points for Plots 1-2 (default 200)
#' @param n         Numeric vector of sample sizes (same length as p).
#'                  Needed for the diagnostic U computation in Plot 3.
#'                  If NULL, Plot 3 is skipped.
#' @param q         Numeric vector of parameter dimensions (same length as p).
#'                  Needed for the diagnostic U computation in Plot 3.
#'                  If NULL, Plot 3 is skipped.
#' @param ...       Additional arguments passed to plot()
#' @return Invisibly, list with:
#'   alpha_grid     — alpha values used in the grid
#'   Cstar_alpha    — C*(alpha) at each grid point
#'   proportions    — achieved proportion at each grid point
#'   Cstar_at_alpha — C* at the pre-chosen alpha (used for Plot 3)
calibration_plot_new <- function(p, ejab, up = 0.1, alpha = 0.05,
                                  grid_range = c(0, 1), grid_n = 200,
                                  n_alpha = 200, n = NULL, q = NULL, ...) {

  # --- Step 1: Estimate C*(alpha) for each alpha on the grid ---
  # This produces the data for Plots 1 and 2
  fit <- estimate_Cstar_alpha(p, ejab, up, grid_range, grid_n, n_alpha)

  # Look up C*(alpha) at the pre-chosen alpha by finding the nearest grid point
  nearest <- which.min(abs(fit$alpha_grid - alpha))
  Cstar_at_alpha <- fit$Cstar_alpha[nearest]

  # --- Plot 1: Calibration curve (Figure 3a style) ---
  # x = alpha grid, y = observed proportion using adaptive C*(alpha)
  # Reference line: slope = 1/up (the theoretical T1E rate)
  # If well-calibrated, the curve tracks the reference line closely
  plot(fit$alpha_grid, fit$proportions,
       type = "l", lwd = 2,
       xlab = expression(alpha),
       ylab = "Observed Proportion",
       xlim = c(0, up),
       ylim = c(0, max(fit$proportions, 1)),
       main = expression(paste("Calibration using ", C^"*", "(", alpha, ")")),
       ...)
  abline(0, 1 / up, col = "red", lty = 2, lwd = 2)

  # --- Plot 2: C*(alpha) vs alpha ---
  # Shows how the optimal threshold varies with alpha
  # Horizontal grey line at C=1 (the natural Bayes factor boundary)
  plot(fit$alpha_grid, fit$Cstar_alpha,
       type = "l", lwd = 2,
       xlab = expression(alpha),
       ylab = expression(C^"*" * (alpha)),
       main = expression(paste(C^"*", "(", alpha, ") vs ", alpha)),
       ...)
  abline(h = 1, col = "grey50", lty = 3)

  # --- Plot 3: Diagnostic QQ-plot at the pre-chosen alpha ---
  # Detect contradictions at the pre-chosen alpha using C*(alpha) from Plot 2,
  # compute the diagnostic U_i values, and produce a QQ-plot.
  # Requires n and q vectors to compute diagnostic_U().
  if (!is.null(n) && !is.null(q)) {
    # Find indices of contradictions: p < alpha AND ejab > C*(alpha)
    idx <- which(p < alpha & ejab > Cstar_at_alpha)
    if (length(idx) > 0) {
      # Compute diagnostic U_i for each contradiction
      U <- diagnostic_U(p[idx], n[idx], q[idx], alpha, Cstar_at_alpha)
      # QQ-plot with best-fit line and simultaneous confidence band
      diagnostic_qqplot_fit(U, alpha = alpha, Cstar = Cstar_at_alpha, ...)
    } else {
      message("No candidates at alpha = ", alpha, " with C* = ",
              round(Cstar_at_alpha, 4))
    }
  }

  invisible(list(alpha_grid = fit$alpha_grid,
                 Cstar_alpha = fit$Cstar_alpha,
                 proportions = fit$proportions,
                 Cstar_at_alpha = Cstar_at_alpha))
}
