f_black_scholes <- function(S, K, T, rf_structure, sigma, 
                            trading_day_convention = 250, 
                            rate_day_convention = 360) {
  ## Black-Scholes European option pricing formula
  ## 
  ## Parameters:
  ##   S (numeric): Current spot price of underlying asset
  ##   K (numeric): Strike price of option
  ##   T (numeric): Time to maturity in years (trading-day convention)
  ##   rf_structure (data.frame): Term structure with columns y_maturity and rate
  ##   sigma (numeric): Implied volatility (annualized)
  ##   trading_day_convention (numeric): Days per year for option convention (default 250)
  ##   rate_day_convention (numeric): Days per year for rate convention (default 360)
  ##
  ## Returns:
  ##   data.frame: 2-row frame with Call and Put prices plus parameters
  
  # Convert T from trading-day years to rate-convention years for interpolation
  T_days <- T * trading_day_convention
  T_rate <- T    
  
  r <- f_interpolate_rates(rf_structure, T_rate)              # interpolate 
  
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

f_interpolate_rates <- function(rf_structure, T) {
  ## Linear interpolation of risk-free rate from term structure
  ##
  ## Parameters:
  ##   rf_structure (data.frame): Term structure with columns y_maturity (years) and rate
  ##   T (numeric): Maturity for interpolation (in years)
  ##
  ## Returns:
  ##   numeric: Interpolated risk-free rate; extrapolates if T outside range
  ##
  ## Logic:
  ##   - If T <= min maturity: returns minimum maturity rate
  ##   - If T >= max maturity: returns maximum maturity rate
  ##   - Otherwise: linear interpolation between surrounding maturities
  if (T <= min(rf_structure$y_maturity)) {
    return(rf_structure$rate[which.min(rf_structure$y_maturity)])
  } else if (T >= max(rf_structure$y_maturity)) {
    return(rf_structure$rate[which.max(rf_structure$y_maturity)])
  } else {
    lower_index <- max(which(rf_structure$y_maturity <= T))
    upper_index <- min(which(rf_structure$y_maturity >= T))
    
    # If T exactly matches a maturity point, return that rate
    if (lower_index == upper_index) {
      return(rf_structure$rate[lower_index])
    }
    
    r_lower <- rf_structure$rate[lower_index]
    r_upper <- rf_structure$rate[upper_index]
    T_lower <- rf_structure$y_maturity[lower_index]
    T_upper <- rf_structure$y_maturity[upper_index]
    
    # Linear interpolation
    r_interpolated <- r_lower + (r_upper - r_lower) * (T - T_lower) / (T_upper - T_lower)
    
    return(r_interpolated)
  }
}
