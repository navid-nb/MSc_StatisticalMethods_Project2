

f_black_scholes <- function(S, K, T, rf_structure, sigma, day_convention = 250) {
## S: Current stock price
## K: Strike price
## T: Time to maturity (in years)
## r: Risk-free interest rate
## sigma: Volatility of the underlying asset
## day_convention: Number of trading days in a year (default is 250)
##rf_structure: Data frame containing the risk-free rate structure with columns 'y_maturity' and 'rate'
  
  r = f_interpolate_rates(rf_structure, T)
  
  d1 <- (log(S / K) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  #Caluclate call price using Black-Scholes formula
  price_call <- S * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  
  #Caclulate put price using put-call parity 
  price_put <- price_call - S + K * exp(-r * T)
  
  #print(paste("Call option price (T =", T * day_convention, ", K =", K, "):", round(price_call, 3), "$"))
  #print(paste("Put option price (T =", T * day_convention, ", K =", K, "):", round(price_put, 3), "$"))
  
  results <- data.frame(
    Option_Type = c("Call", "Put"),
    Maturity_days = c(T, T) * day_convention,
    Strike_Price = c(K, K),
    Risk_Free_Rate = c(r, r),
    Volatility = c(sigma, sigma),
    Price = c(price_call, price_put)
  )
  
  return(results)}

f_interpolate_rates <- function(rf_structure, T) {
  # Interpolate the risk-free rate for the given maturity T
  if (T <= min(rf_structure$y_maturity)) {
    return(rf_structure$rate[which.min(rf_structure$y_maturity)])
  } else if (T >= max(rf_structure$y_maturity)) {
    return(rf_structure$rate[which.max(rf_structure$y_maturity)])
  } else {
    lower_index <- max(which(rf_structure$y_maturity <= T))
    upper_index <- min(which(rf_structure$y_maturity >= T))
    
    r_lower <- rf_structure$rate[lower_index]
    r_upper <- rf_structure$rate[upper_index]
    T_lower <- rf_structure$y_maturity[lower_index]
    T_upper <- rf_structure$y_maturity[upper_index]
    
    # Linear interpolation
    r_interpolated <- r_lower + (r_upper - r_lower) * (T - T_lower) / (T_upper - T_lower)
    
    return(r_interpolated)
  }
}
