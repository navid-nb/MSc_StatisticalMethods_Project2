library(zoo)
library(forecast)
library(quantmod)
library(stats)
library(dplyr)
library(tidyr)
library(PerformanceAnalytics)
library(cvar)
library(MASS)

#Black-Scholes formula


BS_call <- function(K, S, rf, sigma, T_to_mat) {
  #Objective: Compute the price of a call through Black-Scholes
  
  #Parameters:
  
  d1 <- (log(S/K) + T_to_mat * (rf + (sigma^2 /2))) / (sigma * sqrt(T_to_mat))
  d2 <- d1 - (sigma * sqrt(T_to_mat))
  
  call_price <- S * pnorm(d1) - K * exp(-rf * T_to_mat) * pnorm(d2)
  return(call_price)
}

#Volatility Surface setup

vol_surface <- function(alpha, m, tau) {
# Objective: Setup function of implied volatility (parametric surface given)
  
# Parameters:
## alpha = Vector of initial guesses for all 4 alphas in the parametric surface
## m = Moneyness
## tau = Time to maturity
  alpha[1] + 
    alpha[2] * (m - 1)^2 + 
    alpha[3] * (m - 1)^3 + 
    alpha[4] * sqrt(tau)
}


min_dist_vol_surf <- function(alpha) {
# Objective: Setup of objective function to be minimized; the distance between IV of the market and the model.

# Parameters:
## alpha = Vector of initial guesses for all 4 alphas in the parametric surface  

  IV_model <- vol_surface(alpha, m, tau)
  sum(abs(IV_market - IV_model))
}
