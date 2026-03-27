f_vol_surface <- function(alpha, m, tau) {
#alpha is a vector of parameters to be estimated, m is the moneyness (K/S), 
#and tau is the time to maturity in years.
  alpha[1] + 
    alpha[2] * (m - 1)^2 + 
    alpha[3] * (m - 1)^3 + 
    alpha[4] * sqrt(tau)
}

f_objective <- function(alpha, data) {
#Minimize the sum of absolute differences between the model-implied volatilities
# and the observed implied volatilities in the data.
  sigma_model <- f_vol_surface(
    alpha = alpha,
    m     = data$Moneyness,
    tau   = data$Tau
  )
  sum(abs(data$Implied_Volatility - sigma_model))
}