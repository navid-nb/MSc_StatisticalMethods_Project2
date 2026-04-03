f_normal_scenarios_bivar <- function(mu_returns,
                               sigma_returns,
                               n_scenarios,
                               T,
                               spot_price,
                               corr_matrix = NULL) {
  ## Multivariate Gaussian scenarios for levels at horizon T
  ##
  ## Assumption: daily invariants (log-returns) are iid multivariate Gaussian
  ##   R_t ~ N(mu, Sigma_daily)
  ## so the T-day cumulative log-return is
  ##   sum_{t=1..T} R_t ~ N(T * mu, T * Sigma_daily)
  ##
  ## Inputs must be vectors of the same length d (number of drivers).
  ## Returns an (n_scenarios x d) matrix of simulated levels at horizon.

  mu_returns <- as.numeric(mu_returns)
  sigma_returns <- as.numeric(sigma_returns)
  spot_price <- as.numeric(spot_price)

  d <- length(mu_returns)
  stopifnot(length(sigma_returns) == d)
  stopifnot(length(spot_price) == d)

  if (is.null(corr_matrix)) {
    corr_matrix <- diag(1, d)
  } else {
    corr_matrix <- as.matrix(corr_matrix)
    stopifnot(nrow(corr_matrix) == d, ncol(corr_matrix) == d)
  }

  Sigma_daily <- diag(sigma_returns) %*% corr_matrix %*% diag(sigma_returns)
  Sigma_T <- Sigma_daily * T
  mu_T <- mu_returns * T

  chol_Sigma_T <- chol(Sigma_T + diag(1e-12, d))
  z <- matrix(rnorm(n_scenarios * d), nrow = n_scenarios, ncol = d)
  cum_log_returns <- sweep(z %*% chol_Sigma_T, 2, mu_T, "+")

  sweep(exp(cum_log_returns), 2, spot_price, "*")
}

