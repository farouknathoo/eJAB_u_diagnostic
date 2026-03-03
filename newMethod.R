
# EJAB01

ejab01 <- function(p, n, q) {
  stopifnot(p > 0, p < 1, n > 0, q >= 1)

  chi_quant <- qchisq(1 - p, df = q)
  sqrt(n) * exp(-0.5 * (n^(1/q) - 1) * chi_quant)
}




# Empirical Prob, Phat
phat <- function(p, ejab, alpha, C, up) {
  idx <- p <= up
  p_sub <- p[idx]
  ejab_sub <- ejab[idx]

  mean(p_sub <= alpha & ejab_sub > C)
}





# Objective fn for C
objective_C <- function(C, p, ejab, up) {
  p_sub <- p[p <= up]
  p_grid <- sort(unique(p_sub))
  N <- length(p_sub)

  obj <- 0

  for (j in seq_along(p_grid)) {
    alpha <- p_grid[j]
    ph <- phat(p, ejab, alpha, C, up)
    obj <- obj + (ph - alpha / up)^2
  }

  obj
}





# Grid search for C
estimate_Cstar <- function(p, ejab, up, grid = seq(1/3, 3, length.out = 200)) {
  vals <- sapply(grid, objective_C, p = p, ejab = ejab, up = up)

  idx <- which.min(vals)
  list(
    Cstar = grid[idx],
    objective = vals[idx]
  )
}





# Type 1 identifier
detect_type1 <- function(p, ejab, alpha, Cstar) {
  which(p <= alpha & ejab > Cstar)
}





# Ui
diagnostic_U <- function(p, n, q, alpha, Cstar) {
  denom_term <- function(ni, qi) {
    1 - pchisq(
      (2 * ni^(1/qi) / (ni^(1/qi) - 1)) * log(sqrt(ni)/Cstar),
      df = qi
    )
  }

  sapply(seq_along(p), function(i) {
    d <- denom_term(n[i], q[i])
    (p[i] - d) / (alpha - d)
  })
}





# Wrapper
ejab_pipeline <- function(df, up, alpha) {
  ejab <- ejab01(df$p, df$N, df$q)
  Cfit <- estimate_Cstar(df$p, ejab, up)

  idx <- detect_type1(df$p, ejab, alpha, Cfit$Cstar)

  list(
    Cstar = Cfit$Cstar,
    objective = Cfit$objective,
    candidates = df[idx, ],
    ejab = ejab
  )
}








