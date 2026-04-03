f_normal_scenarios <- function(mu_returns, sigma_returns, n_scenarios, T, spot_price) {
  ## Monte Carlo scenario generation under univariate Gaussian model (GBM)
  ##
  ## Parameters:
  ##   mu_returns (numeric): Mean daily log-return
  ##   sigma_returns (numeric): Standard deviation of daily log-return
  ##   n_scenarios (numeric): Number of Monte Carlo scenarios to generate
  ##   T (numeric): Number of days ahead for projection
  ##   spot_price (numeric): Current price of underlying asset
  ##
  ## Returns:
  ##   numeric vector: Simulated price scenarios at time T
  
  # Convert spot_price to numeric (handles xts/array inputs)
  spot_price <- as.numeric(spot_price)
  
  # Generate random log-returns for T days
  log_return_scenarios <- matrix(rnorm(T * n_scenarios, mean = mu_returns, sd = sigma_returns), nrow = n_scenarios, ncol = T)
  
  # Calculate the cumulative log-returns for each scenario
  cumulative_log_returns <- rowSums(log_return_scenarios)
  
  # Calculate the price scenarios by "exponentiating" the cumulative log-returns and multiplying by the spot price
  price_scenarios <- spot_price * exp(cumulative_log_returns)
  
  return(price_scenarios)
}
