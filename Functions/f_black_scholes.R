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
  ## Linear interpolation of risk-free rate from term structure (handles scalar or vector T)
  ##
  ## Parameters:
  ##   rf_structure (data.frame): Term structure with columns y_maturity (years) and rate
  ##   T (numeric): Maturity for interpolation (scalar or vector, in years)
  ##
  ## Returns:
  ##   numeric: Interpolated risk-free rate(s); extrapolates if T outside range
  ##
  ## Logic:
  ##   - Vectorized: works with scalar T or vector of T values
  ##   - If T <= min maturity: returns minimum maturity rate
  ##   - If T >= max maturity: returns maximum maturity rate
  ##   - Otherwise: linear interpolation between surrounding maturities
  
  # Handle vector input
  if (length(T) > 1) {
    return(sapply(T, function(t) f_interpolate_rates(rf_structure, t)))
  }
  
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

f_black_scholes_vectorized <- function(S, K, T, rf_structure, sigma,
                                       trading_day_convention = 250,
                                       rate_day_convention = 360) {
  ## Vectorized Black-Scholes for single spot price, multiple strikes/maturities
  ## Optimized for parallel Monte Carlo: returns call prices as vector
  ##
  ## Parameters:
  ##   S (numeric): Single spot price
  ##   K (numeric vector): Strike prices
  ##   T (numeric vector): Times to maturity in years
  ##   rf_structure (data.frame): Term structure
  ##   sigma (numeric): Implied volatility (scalar or vector, length = length(K))
  ##   trading_day_convention (numeric): Default 250
  ##   rate_day_convention (numeric): Default 360
  ##
  ## Returns:
  ##   numeric vector: Call prices (same length as K)

  # Ensure vectors
  K <- as.vector(K)
  T <- as.vector(T)
  sigma <- as.vector(sigma)
  
  # Vectorized time convention conversion
  T_days <- T * trading_day_convention
  T_rate <- T
  
  # Vectorized rate interpolation
  r <- f_interpolate_rates(rf_structure, T_rate)
  
  # Vectorized Black-Scholes calculations
  d1 <- (log(S / K) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  price_call <- S * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  
  return(price_call)
}

f_calculate_forward_rates <- function(rf_structure, days_passed, trading_days_in_year) {
  ## Convert t=0 spot rate term structure to forward rates from (t=0 + days_passed) onward
  ##
  ## Parameters:
  ##   rf_structure (data.frame): Term structure at t=0 with columns y_maturity, rate
  ##   days_passed (numeric): Days elapsed (e.g., 5) to forward from
  ##   trading_days_in_year (numeric): Day convention (typically 250)
  ##
  ## Returns:
  ##   data.frame: Forward rate term structure with new y_maturity (remaining time) and forward rates
  ##
  ## Formula (continuously compounded):
  ##   f(t,T) = [r(T) * T - r(t) * t] / (T - t)
  
  t <- days_passed / trading_days_in_year
  r_t <- f_interpolate_rates(rf_structure, t)
  
  # Forward rates are only defined for maturities strictly greater than the horizon t.
  idx <- rf_structure$y_maturity > t
  rf_long <- rf_structure[idx, , drop = FALSE]
  
  forward_rates <- (rf_long$rate * rf_long$y_maturity - r_t * t) /
                   (rf_long$y_maturity - t)
  
  data.frame(
    y_maturity = rf_long$y_maturity - t,
    rate = forward_rates
  )
}
