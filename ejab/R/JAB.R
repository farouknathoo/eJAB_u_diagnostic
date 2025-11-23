#' JAB01: Compute the Jeffreys Approximate Objective Bayes factor from p-value
#'
#' This function calculates the Jeffreys Approximate Objective Bayes (JAB) factor
#' based on the formula:
#' N^(k/2) * exp(-0.5 * (N - 1)/N * Q_k(1 - p-value)).
#'
#' @param N Sample size (number of observations).
#' @param pVal P-value.
#' @param model The statistical model being used, one of:
#'   - "t-test"
#'   - "linear_regression"
#'   - "logistic_regression"
#'   - "cox"
#'   - "anova" (one-way ANOVA)
#'   - "repeated_measures"
#'   - "chi_squared"
#' @param R (Optional) Number of rows for chi-squared tests.
#' @param C (Optional) Number of columns for chi-squared tests.
#' @param I (Optional) Number of conditions for one-way ANOVA and repeated measures.
#' @return A numeric value corresponding to the Jeffreys Approximate Bayes factor.
#' @examples
#' JAB01(100, 0.05, "t-test")
#' JAB01(200, 0.01, "chi_squared", R = 3, C = 4)
#' @export
JAB01 <- function(N, pVal, model, R = NULL, C = NULL, I = NULL) {
  # Validate inputs
  if (!is.numeric(N) || N <= 0) stop("N must be a positive numeric value.")
  if (!is.numeric(pVal) || pVal <= 0 || pVal >= 1) stop("pVal must be between 0 and 1.")
  if (!model %in% c("t-test", "linear_regression", "logistic_regression", "cox",
                    "anova", "repeated_measures", "chi_squared","kruskal_wallis","wilcoxon","mann_whitney")) {
    stop("Invalid model type. Choose from 't-test', 'linear_regression', 'logistic_regression', 'cox', 'anova', 'repeated_measures', 'chi_squared','kruskal_wallis','mann_whitney'.")
  }

  # Determine k based on the model type
  k <- switch(
    model,
    "t-test" = 1,
    "wilcoxon" = 1,
    "mann_whitney" = 1,
    "linear_regression" = 1,
    "logistic_regression" = 1,
    "cox" = 1,
    "anova" = if (!is.null(I)) I - 1 else stop("I must be specified for 'anova'."),
    "kruskal_wallis" = if (!is.null(I)) I - 1 else stop("I must be specified for 'kruskal_wallis'."),
    "repeated_measures" = if (!is.null(I)) I - 1 else stop("I must be specified for 'repeated_measures'."),
    "chi_squared" = if (!is.null(R) && !is.null(C)) (R - 1) * (C - 1) else stop("R and C must be specified for 'chi_squared'.")
  )

  # Compute JAB01
  #JAB <- N^(k / 2) * exp(-0.5 * ((N - 1) / N) * stats::qchisq(1 - pVal, df = k))
  JAB <- sqrt(N) * exp(-0.5 * (((N^(1/k) - 1) / N^(1/k)) * stats::qchisq(1 - pVal, df = k)))

  return(JAB)
}
