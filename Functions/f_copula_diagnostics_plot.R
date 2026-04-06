## Copula Fit Diagnostics Plots
##
## Visualizes the copula-marginal model fitting process and results.
## Includes fitted Student-t distributions, joint/marginal uniforms,
## and original returns for model validation.
##
## Args:
##   daily_log_returns: numeric vector of asset log-returns
##   daily_log_returns_vix: numeric vector of VIX log-returns
##   df_asset: degrees of freedom for asset Student-t
##   df_vix: degrees of freedom for VIX Student-t
##   mu_asset: estimated location parameter for asset
##   mu_vix: estimated location parameter for VIX
##   scale_asset: estimated scale parameter for asset
##   scale_vix: estimated scale parameter for VIX
##   sigma_asset: standard deviation for asset (derived from scale)
##   sigma_vix: standard deviation for VIX (derived from scale)
##   uniforms: data frame with U_Asset and U_VIX columns
##   copula_correlation: estimated Gaussian copula correlation
##   output_path: path to save PNG output (optional)
##
## Returns:
##   ggplot object with combined diagnostic plots

plot_copula_fit_diagnostics <- function(daily_log_returns, 
                                        daily_log_returns_vix,
                                        df_asset, 
                                        df_vix,
                                        mu_asset, 
                                        mu_vix,
                                        scale_asset, 
                                        scale_vix,
                                        sigma_asset, 
                                        sigma_vix,
                                        uniforms,
                                        copula_correlation,
                                        output_path = NULL) {
  
  # Plot 1: Asset returns with fitted Student-t distribution overlaid
  p_asset_fit <- ggplot(data.frame(returns = daily_log_returns), aes(x = returns)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#3498DB", alpha = 0.6, color = NA) +
    stat_function(fun = function(x) dt((x - mu_asset) / scale_asset, df = df_asset) / scale_asset,
                  color = "#E74C3C", linewidth = 1, linetype = "solid") +
    labs(title = "Asset Returns with Fitted Student-t(df=10)",
         x = "Log-returns", y = "Density") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11))
  
  # Plot 2: VIX returns with fitted Student-t distribution overlaid
  p_vix_fit <- ggplot(data.frame(returns = daily_log_returns_vix), aes(x = returns)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#E74C3C", alpha = 0.6, color = NA) +
    stat_function(fun = function(x) dt((x - mu_vix) / scale_vix, df = df_vix) / scale_vix,
                  color = "#2C3E50", linewidth = 1, linetype = "solid") +
    labs(title = "VIX Returns with Fitted Student-t(df=5)",
         x = "Log-returns", y = "Density") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11))
  
  # Plot 3: Marginal distribution of Asset uniforms
  p1 <- ggplot(uniforms, aes(x = U_Asset)) +
    geom_histogram(bins = 50, fill = "#3498DB", alpha = 0.7, color = NA) +
    labs(title = "Asset Marginal (Uniform)", x = "U_Asset", y = "Frequency") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10))
  
  # Plot 4: Marginal distribution of VIX uniforms
  p2 <- ggplot(uniforms, aes(x = U_VIX)) +
    geom_histogram(bins = 50, fill = "#E74C3C", alpha = 0.7, color = NA) +
    labs(title = "VIX Marginal (Uniform)", x = "U_VIX", y = "Frequency") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10))
  
  # Plot 5: Original bivariate returns
  p3 <- ggplot(data.frame(Asset = daily_log_returns, VIX = daily_log_returns_vix),
               aes(x = Asset, y = VIX)) +
    geom_point(color = rgb(0, 0, 1, 0.3), size = 0.8, alpha = 0.6) +
    labs(title = "Original Returns (Asset vs VIX)",
         x = "Asset log-returns", y = "VIX log-returns") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))
  
  # Plot 6: Bivariate uniform samples 
  p4 <- ggplot(uniforms, aes(x = U_Asset, y = U_VIX)) +
    geom_point(color = rgb(1, 0, 0, 0.3), size = 0.8, alpha = 0.6) +
    labs(title = "Pseudo Data Generated (Joint Uniforms)",
         x = "U_Asset", y = "U_VIX") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))
  
  # Combined plot: Fitted distributions at top, then joint/uniforms below
  copula_plot <- (p_asset_fit + p_vix_fit) / (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title    = "Copula-Marginal Model: Student-t Marginals with Gaussian Copula",
      subtitle = paste0("Asset: t(df=", df_asset, ", μ=", round(mu_asset, 4), ", σ=", round(sigma_asset, 4), ") | ",
                        "VIX: t(df=", df_vix, ", μ=", round(mu_vix, 4), ", σ=", round(sigma_vix, 4), ") | ",
                        "Copula ρ=", round(copula_correlation, 4)),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "#566573", size = 9),
        plot.caption  = element_text(color = "#566573", size = 8, hjust = 0)
      )
    )
  
  print(copula_plot)
  
  if (!is.null(output_path)) {
    ggsave(filename = output_path, plot = copula_plot, width = 18, height = 16, dpi = 300)
  }
  
  return(copula_plot)
}
