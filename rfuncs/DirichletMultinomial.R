posterior_analysis <- function(data_summary2, n_sim = 5000, alpha = 0.05) {
  # Prepare the study data
  study.dat <- data_summary2[!is.na(data_summary2$JAB), ]

  # Get unique study identifiers
  pmid <- unique(study.dat$nctId)
  M.S <- length(pmid)

  # Precompute conditions for all studies
  alpha_params <- list()
  study_sizes <- numeric(M.S)

  for (i in 1:M.S) {
    ID <- pmid[i]
    temp.dat <- study.dat[study.dat$nctId == ID, ]
    study_sizes[i] <- nrow(temp.dat)

    # Count outcomes in each category
    y1 <- sum((temp.dat$pValue >= alpha) & (temp.dat$JAB < 1/3))
    y2 <- sum((temp.dat$pValue >= alpha) & (temp.dat$JAB >= 1/3))
    y3 <- sum((temp.dat$pValue < alpha) & (temp.dat$JAB < 1/3))
    y4 <- sum((temp.dat$pValue < alpha) & (temp.dat$JAB >= 1/3))

    # Store Dirichlet parameters
    alpha_params[[i]] <- c(1 + y1, 1 + y2, 1 + y3, 1 + y4)
  }

  # Initialize results vector
  p.posterior <- numeric(n_sim)

  # For each study, draw n_sim samples at once from Dirichlet
  dirichlet_samples <- list()
  for (i in 1:M.S) {
    dirichlet_samples[[i]] <- gtools::rdirichlet(n_sim, alpha_params[[i]])
  }

  # Calculate posterior probabilities for all simulations at once
  for (j in 1:n_sim) {
    s.num <- 0
    s.denom <- 0

    for (i in 1:M.S) {
      s.num <- s.num + study_sizes[i] * dirichlet_samples[[i]][j, 4]
      s.denom <- s.denom + study_sizes[i] * (dirichlet_samples[[i]][j, 3] + dirichlet_samples[[i]][j, 4])
    }

    p.posterior[j] <- if(s.denom == 0) NA else s.num / s.denom
  }

  # Calculate summary statistics
  posterior_mean <- mean(p.posterior, na.rm = TRUE)
  posterior_se <- sqrt(var(p.posterior, na.rm = TRUE))

  return(list(mean = posterior_mean, se = posterior_se, data = p.posterior ))
}
