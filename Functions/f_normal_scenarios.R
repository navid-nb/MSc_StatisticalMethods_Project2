f_normal_scenarios <- function(mu_returns, sigma_returns, n_scenarios, T, spot_price) {
  # Generate iid normally distributed scenarios for the underlying asset price
  # mu_returns: Mean of the daily log-returns
  # sigma_returns: Standard deviation of the daily log-returns
  # n_scenarios: Number of scenarios to generate
  # T: days ahead for which the scenarios are generated 
  # spot_price: Current price of the underlying asset
  
  # Generate random log-returns for T days
  log_return_scenarios <- matrix(rnorm(T * n_scenarios, mean =mu_returns, sd = sigma_returns), nrow = n_scenarios, ncol = T)
  
  # Calculate the cumulative log-returns for each scenario
  cumulative_log_returns <- rowSums(log_return_scenarios)
  
  # Calculate the price scenarios by "exponentiating" the cumulative log-returns and multiplying by the spot price
  price_scenarios <- spot_price * exp(cumulative_log_returns)
  
  return(price_scenarios)
}
