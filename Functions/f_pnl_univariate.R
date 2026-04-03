## Univariate Gaussian Model - P&L Calculation

## Calculate P&L Distribution - Univariate Model
## 
## Computes P&L distribution for portfolio assuming only spot price varies
## (implied volatility held constant). Applies vectorized Black-Scholes pricing
## across 10,000 Monte Carlo scenarios.
##
## Args:
##   scenarios_asset: numeric vector of simulated spot prices at time T
##   portfolio: data.frame with columns Strike_Price, Maturity_days
##   rf_structure: data.frame with columns y_maturity, rate (term structure)
##   trading_days_in_year: day convention for option maturity (default 250)
##   last_volatility: current implied volatility (scalar)
##   portfolio_value_inception: initial portfolio value (scalar)
##
## Returns:
##   list with elements:
##     - pnl_distribution: vector of 10,000 P&L values
##     - VaR95: Value at Risk at 95% confidence level
##     - ES95: Expected Shortfall at 95% confidence level

calculate_pnl_univariate <- function(scenarios_asset, portfolio, rf_structure, 
                                     trading_days_in_year, last_volatility, 
                                     portfolio_value_inception) {
  
  portfolio_value_five_days_ahead <- sapply(scenarios_asset, function(scenario_price) {
    call_prices <- f_black_scholes_vectorized(
      S = scenario_price,
      K = portfolio$Strike_Price,
      T = (portfolio$Maturity_days - 5) / trading_days_in_year,
      rf_structure = rf_structure,
      sigma = rep(last_volatility, nrow(portfolio))
    )
    sum(call_prices)
  })
  
  pnl_distribution <- portfolio_value_five_days_ahead - portfolio_value_inception
  VaR95 <- quantile(pnl_distribution, 0.05)
  ES95  <- mean(pnl_distribution[pnl_distribution <= VaR95])
  
  return(list(
    pnl_distribution = pnl_distribution,
    VaR95 = VaR95,
    ES95 = ES95
  ))
}

## Univariate P&L Histogram with Risk Metrics

## Plot P&L Distribution - Univariate Model
##
## Creates histogram with density overlay and risk metric annotations.
## Displays Value at Risk (95% confidence, cyan line) and Expected Shortfall 
## (95% confidence, dark blue line) on the distribution.
##
## Args:
##   pnl_distribution: numeric vector of P&L values
##   VaR95: scalar, 95% confidence VaR (5th percentile)
##   ES95: scalar, 95% confidence Expected Shortfall
##   output_path: character string for PNG output location (optional, if NULL no file saved)
##
## Returns:
##   ggplot object (invisibly prints to console)

plot_pnl_univariate <- function(pnl_distribution, VaR95, ES95, output_path = NULL) {
  
  pnl_df <- data.frame(PnL = pnl_distribution)
  
  pnl_distrib_plot <- ggplot(pnl_df, aes(x = PnL)) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = 1, fill = "#2C3E50", color = "#AEB6BF", alpha = 0.85) +
    geom_density(color = "#E74C3C", linewidth = 0.8, adjust = 1) +
    geom_vline(xintercept = VaR95, color = "cyan", linewidth = 0.9) +
    geom_vline(xintercept = ES95,  color = "darkblue", linewidth = 0.9) +
    annotate("text", x = VaR95, y = Inf, 
             label = paste0("VaR 95%: $", round(VaR95, 2)),
             color = "cyan", hjust = -0.1, vjust = 2, size = 3.5, fontface = "bold") +
    annotate("text", x = ES95, y = Inf,  
             label = paste0("ES 95%: $",  round(ES95, 2)),
             color = "darkblue", hjust = -0.1, vjust = 4, size = 3.5, fontface = "bold") +
    labs(title    = "P&L Distribution — Portfolio of European Call Options",
         subtitle = paste0("10,000 scenarios | 5-day horizon | ",
                           "VaR 95%: $", round(VaR95, 2), 
                           " | ES 95%: $", round(ES95, 2)),
         x = "P&L ($)", 
         y = "Density",
         caption = "Univariate Gaussian model — one risk driver | Red curve: KDE (Silverman bandwidth)") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "#566573", size = 10),
      plot.caption  = element_text(color = "#566573", size = 9, hjust = 0),
      panel.grid.minor = element_blank()
    )
  
  print(pnl_distrib_plot)
  
  if (!is.null(output_path)) {
    ggsave(filename = output_path, plot = pnl_distrib_plot, width = 8, height = 5, dpi = 300)
  }
  
  return(pnl_distrib_plot)
}

