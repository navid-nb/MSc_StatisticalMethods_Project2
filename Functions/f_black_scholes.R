f_black_scholes <- function(S, K, T, rf_structure, sigma, 
                            trading_day_convention = 250, 
                            rate_day_convention = 360) {
  ## T: Time to maturity in years (using trading_day_convention)
  ## trading_day_convention: For option maturity (default 250)
  ## rate_day_convention: For term structure interpolation (default 365)
  
  # Convert T from trading-day years to rate-convention years for interpolation
  T_days <- T * trading_day_convention                        # back to raw days
  T_rate <- T_days / rate_day_convention                      # into rate convention
  
  r <- f_interpolate_rates(rf_structure, T_rate)              # interpolate in rate convention
  
  # Black-Scholes uses T in trading-day convention throughout
  d1 <- (log(S / K) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  price_call <- S * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  price_put  <- price_call - S + K * exp(-r * T)
  
  results <- data.frame(
    Option_Type    = c("Call", "Put"),
    Maturity_days  = c(T_days, T_days),
    Strike_Price   = c(K, K),
    Risk_Free_Rate = c(r, r),
    Volatility     = c(sigma, sigma),
    Price          = c(price_call, price_put)
  )
  
  return(results)
}

f_interpolate_rates <- function(rf_structure, T, day_convention) {
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
