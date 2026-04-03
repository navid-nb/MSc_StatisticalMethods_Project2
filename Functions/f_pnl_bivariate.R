## Bivariate Gaussian Model - P&L Calculation
##
## Computes P&L distribution for a portfolio of European call options when
## two risk drivers are simulated jointly:
##   - Underlying spot (S)
##   - Implied volatility proxy (e.g., VIX)
##
## Pricing uses the vectorized Black-Scholes implementation from f_black_scholes.R
## to evaluate all options in the portfolio per Monte Carlo scenario.

calculate_pnl_bivariate <- function(scenarios_asset,
                                    scenarios_vix,
                                    portfolio,
                                    rf_structure,
                                    trading_days_in_year,
                                    portfolio_value_inception,
                                    horizon_days = 5) {
  stopifnot(length(scenarios_asset) == length(scenarios_vix))

  # Remaining time-to-maturity at horizon (in years, trading-day convention)
  T_remaining <- (portfolio$Maturity_days - horizon_days) / trading_days_in_year

  portfolio_value_horizon <- mapply(function(spot_h, vol_h) {
    matured <- T_remaining <= 0
    call_prices <- numeric(nrow(portfolio))

    if (any(matured)) {
      call_prices[matured] <- pmax(spot_h - portfolio$Strike_Price[matured], 0)
    }

    alive <- !matured
    if (any(alive)) {
      call_prices[alive] <- f_black_scholes_vectorized(
        S = spot_h,
        K = portfolio$Strike_Price[alive],
        T = T_remaining[alive],
        rf_structure = rf_structure_forward,
        sigma = rep(vol_h, sum(alive)),
        trading_day_convention = trading_days_in_year
      )
    }

    sum(call_prices)
  }, scenarios_asset, scenarios_vix)

  pnl_distribution <- portfolio_value_horizon - portfolio_value_inception
  VaR95 <- quantile(pnl_distribution, 0.05)
  ES95  <- mean(pnl_distribution[pnl_distribution <= VaR95])

  list(
    pnl_distribution = pnl_distribution,
    VaR95 = VaR95,
    ES95 = ES95
  )
}

plot_pnl_bivariate <- function(pnl_distribution, VaR95, ES95, output_path = NULL) {
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
         caption = "Bivariate Gaussian model — two risk drivers | Red curve: KDE (Silverman bandwidth)") +
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

  pnl_distrib_plot
}
