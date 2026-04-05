## Volatility surface workflow helpers
##
## Assumptions in this repo (kept consistent with existing code):
## - Option maturities for pricing use trading_days_in_year (typically 250)
## - The term structure passed to BMS is whatever the caller provides (e.g., spot curve at T,
##   forward curve at T+5 as built elsewhere)
## - Adjusted volatility follows the assignment formula:
##     sigma_tilde(m, tau) = sigma_hat(m, tau) - (alpha1 + alpha4) + VIX_level

f_vol_surface <- function(alpha, m, tau) {
  # Parametric volatility surface (assignment formula)
  # alpha: vector (alpha1..alpha4)
  # m: moneyness K/S
  # tau: time-to-maturity in years
  alpha[1] +
    alpha[2] * (m - 1)^2 +
    alpha[3] * (m - 1)^3 +
    alpha[4] * sqrt(tau)
}

f_objective <- function(alpha, data) {
  # L1 calibration objective: minimize absolute distance to market IV quotes
  sigma_model <- f_vol_surface(
    alpha = alpha,
    m = data$Moneyness,
    tau = data$Tau
  )
  sum(abs(data$Implied_Volatility - sigma_model))
}

prepare_options_data_otm <- function(calls, puts, spot_price) {
  # Build calibration dataset from market calls/puts, then keep only OTM quotes:
  # - OTM calls: K >= S  (moneyness >= 1)
  # - OTM puts : K <= S  (moneyness <= 1)
  spot_price <- as.numeric(spot_price)

  options_data <- dplyr::bind_rows(
    as.data.frame(calls) %>% dplyr::mutate(Type = "Call"),
    as.data.frame(puts)  %>% dplyr::mutate(Type = "Put")
  ) %>%
    dplyr::mutate(Moneyness = K / spot_price) %>%
    dplyr::select(
      Type,
      Strike_Price = K,
      Implied_Volatility = IV,
      Moneyness,
      Tau = tau
    ) %>%
    dplyr::filter(
      (Type == "Call" & Moneyness >= 1) |
        (Type == "Put" & Moneyness <= 1)
    )

  options_data
}

fit_vol_surface_otm <- function(options_data, last_volatility) {
  # Fit alpha by minimizing absolute distance to market IV (L1 loss).
  a1 <- as.numeric(last_volatility)
  alpha_init <- c(a1, 0.1, 0.1, 0.1)

  optimization <- optim(
    par    = alpha_init,
    fn     = f_objective,
    data   = options_data,
    method = "Nelder-Mead"
  )

  alpha_hat <- optimization$par
  names(alpha_hat) <- c("alpha1", "alpha2", "alpha3", "alpha4")

  list(alpha_hat = alpha_hat, optimization = optimization)
}

plot_vol_surface_3d <- function(options_data, alpha_hat, output_path = NULL) {
  # Simple 3D plot: fitted surface + market IV points.
  #
  # Returns a plotly object; optionally saves to HTML if htmlwidgets is available.
  m_grid <- seq(min(options_data$Moneyness), max(options_data$Moneyness), length.out = 60)
  tau_grid <- seq(min(options_data$Tau), max(options_data$Tau), length.out = 60)

  z <- outer(m_grid, tau_grid, function(m, tau) f_vol_surface(alpha_hat, m, tau))

  p <- plotly::plot_ly() %>%
    plotly::add_surface(
      x = m_grid,
      y = tau_grid,
      z = t(z),
      opacity = 0.65,
      showscale = FALSE,
      name = "Fitted surface"
    ) %>%
    plotly::add_trace(
      data = options_data,
      x = ~Moneyness,
      y = ~Tau,
      z = ~Implied_Volatility,
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = 2,
        color = ~Implied_Volatility,
        colorscale = "Viridis",
        showscale = TRUE,
        colorbar = list(title = "Market IV")
      ),
      name = "Market IV"
    ) %>%
    plotly::layout(
      title = "Fitted IV surface (OTM quotes)",
      scene = list(
        xaxis = list(title = "Moneyness (K/S)"),
        yaxis = list(title = "Tau (years)"),
        zaxis = list(title = "Implied Volatility")
      )
    )

  if (!is.null(output_path) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    htmlwidgets::saveWidget(p, file = output_path, selfcontained = TRUE)
  }

  p
}

f_adjusted_sigma <- function(alpha_hat, moneyness, tau, vix_level) {
  # Assignment adjustment:
  #   sigma_tilde = sigma_hat - (alpha1 + alpha4) + VIX_level
  fitted <- f_vol_surface(alpha_hat, moneyness, tau)
  fitted - (alpha_hat[1] + alpha_hat[4]) + vix_level
}

portfolio_surface_inputs <- function(portfolio,
                                     alpha_hat,
                                     spot_price,
                                     vix_level,
                                     trading_days_in_year,
                                     horizon_days = 0) {
  tau <- (portfolio$Maturity_days - horizon_days) / trading_days_in_year
  m <- portfolio$Strike_Price / spot_price
  fitted <- f_vol_surface(alpha_hat, m, tau)
  adjusted <- fitted - (alpha_hat[1] + alpha_hat[4]) + vix_level

  data.frame(
    Strike_Price = portfolio$Strike_Price,
    Maturity_days = portfolio$Maturity_days,
    Tau = tau,
    Moneyness = m,
    Fitted_IV = fitted,
    Adjusted_Sigma = adjusted
  )
}

portfolio_surface_prices_sigma_hat <- function(portfolio,
                                               alpha_hat,
                                               spot_price,
                                               rf_structure,
                                               trading_days_in_year,
                                               horizon_days = 0) {
  # Reprice the portfolio using the fitted parametric surface sigma_hat(m, tau) directly
  # (i.e., without the VIX level adjustment).
  spot_price <- as.numeric(spot_price)

  strikes <- portfolio$Strike_Price
  tau <- (portfolio$Maturity_days - horizon_days) / trading_days_in_year
  moneyness <- strikes / spot_price

  fitted_iv <- f_vol_surface(alpha_hat, moneyness, tau)

  matured <- tau <= 0
  call_prices <- numeric(length(strikes))

  if (any(matured)) {
    call_prices[matured] <- pmax(spot_price - strikes[matured], 0)
  }

  alive <- !matured
  if (any(alive)) {
    call_prices[alive] <- f_black_scholes_vectorized(
      S = spot_price,
      K = strikes[alive],
      T = tau[alive],
      rf_structure = rf_structure,
      sigma = fitted_iv[alive],
      trading_day_convention = trading_days_in_year
    )
  }

  list(
    option_table = data.frame(
      Strike_Price = strikes,
      Maturity_days = portfolio$Maturity_days,
      Tau = tau,
      Moneyness = moneyness,
      Sigma_hat = fitted_iv,
      Call_Price = call_prices
    ),
    portfolio_value = sum(call_prices)
  )
}

portfolio_value_surface <- function(spot_price,
                                    vix_level,
                                    portfolio,
                                    alpha_hat,
                                    rf_structure,
                                    trading_days_in_year,
                                    horizon_days = 0) {
  strikes <- portfolio$Strike_Price
  tau <- (portfolio$Maturity_days - horizon_days) / trading_days_in_year

  matured <- tau <= 0
  call_prices <- numeric(length(strikes))

  if (any(matured)) {
    call_prices[matured] <- pmax(spot_price - strikes[matured], 0)
  }

  alive <- !matured
  if (any(alive)) {
    m <- strikes[alive] / spot_price
    sigma_adj <- f_adjusted_sigma(alpha_hat, m, tau[alive], vix_level)
    call_prices[alive] <- f_black_scholes_vectorized(
      S = spot_price,
      K = strikes[alive],
      T = tau[alive],
      rf_structure = rf_structure,
      sigma = sigma_adj,
      trading_day_convention = trading_days_in_year
    )
  }

  sum(call_prices)
}

calculate_pnl_vol_surface <- function(portfolio,
                                      alpha_hat,
                                      spot_T,
                                      vix_T,
                                      rf_structure_T,
                                      rf_structure_T_plus_h,
                                      trading_days_in_year,
                                      horizon_days,
                                      spot_scenarios_T_plus_h,
                                      vix_scenarios_T_plus_h) {
  stopifnot(length(spot_scenarios_T_plus_h) == length(vix_scenarios_T_plus_h))

  value_T <- portfolio_value_surface(
    spot_price = spot_T,
    vix_level = vix_T,
    portfolio = portfolio,
    alpha_hat = alpha_hat,
    rf_structure = rf_structure_T,
    trading_days_in_year = trading_days_in_year,
    horizon_days = 0
  )

  value_T_plus_h <- mapply(
    function(spot_h, vix_h) {
      portfolio_value_surface(
        spot_price = spot_h,
        vix_level = vix_h,
        portfolio = portfolio,
        alpha_hat = alpha_hat,
        rf_structure = rf_structure_T_plus_h,
        trading_days_in_year = trading_days_in_year,
        horizon_days = horizon_days
      )
    },
    spot_scenarios_T_plus_h,
    vix_scenarios_T_plus_h
  )

  pnl_distribution <- value_T_plus_h - value_T
  VaR95 <- quantile(pnl_distribution, 0.05)
  ES95  <- mean(pnl_distribution[pnl_distribution <= VaR95])

  list(
    value_T = value_T,
    value_T_plus_h = value_T_plus_h,
    pnl_distribution = pnl_distribution,
    VaR95 = VaR95,
    ES95 = ES95
  )
}

plot_pnl_vol_surface <- function(pnl_distribution,
                                 VaR95,
                                 ES95,
                                 horizon_days,
                                 caption,
                                 output_path = NULL) {
  pnl_df <- data.frame(PnL = pnl_distribution)

  p <- ggplot2::ggplot(pnl_df, ggplot2::aes(x = PnL)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = after_stat(density)),
      binwidth = 0.5, fill = "#27AE60", color = "#AEB6BF", alpha = 0.85
    ) +
    ggplot2::geom_density(color = "#E74C3C", linewidth = 0.8, adjust = 1) +
    ggplot2::geom_vline(xintercept = VaR95, color = "cyan", linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = ES95,  color = "darkblue", linewidth = 0.9) +
    ggplot2::annotate(
      "text", x = VaR95, y = Inf,
      label = paste0("VaR 95%: $", round(VaR95, 2)),
      color = "cyan", hjust = -0.1, vjust = 2, size = 3.5, fontface = "bold"
    ) +
    ggplot2::annotate(
      "text", x = ES95, y = Inf,
      label = paste0("ES 95%: $",  round(ES95, 2)),
      color = "darkblue", hjust = -0.1, vjust = 4, size = 3.5, fontface = "bold"
    ) +
    ggplot2::labs(
      title = "P&L Distribution — Volatility Surface Model",
      subtitle = paste0(
        length(pnl_distribution), " scenarios | ", horizon_days, "-day horizon | ",
        "VaR 95%: $", round(VaR95, 2), " | ES 95%: $", round(ES95, 2)
      ),
      x = "P&L ($)",
      y = "Density",
      caption = caption
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(color = "#566573", size = 10),
      plot.caption = ggplot2::element_text(color = "#566573", size = 9, hjust = 0),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(output_path)) {
    ggplot2::ggsave(filename = output_path, plot = p, width = 8, height = 5, dpi = 300)
  }

  p
}
