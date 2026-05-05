# Statistical Methods — Portfolio P&L Distribution & VaR/ES under Multiple Risk-Factor Models

End-to-end estimation of the 5-day P&L distribution, 95% Value-at-Risk, and Expected Shortfall for a portfolio of four S&P 500 European call options. Four risk-factor models are estimated and benchmarked, each capturing a richer joint structure than the last: univariate Gaussian (spot only), bivariate Gaussian (spot + VIX), copula-marginal model with Student-t marginals and a Gaussian copula, and a parametric volatility-surface model with VIX-shifted dynamics.

---

## Methods & Models

**Portfolio.** Four European calls on the S&P 500 priced via Black-Scholes-Merton: strikes 1600/1650 (20-day) and 1750/1800 (40-day). Portfolio repriced at the 5-day horizon using forward rates derived from the t=0 zero-coupon term structure to ensure correct discounting at T+5.

**1. One Risk Driver — Univariate Gaussian Model**
Closed-form MLE on daily SPX log-returns: $\hat\mu = \bar R$, $\hat\sigma^2 = \tfrac{1}{n}\sum (R_i-\bar R)^2$ (n-divisor for MLE). Time-scaled to a 5-day horizon as $\sum_{t=1}^T R_t \sim \mathcal N(T\hat\mu, T\hat\sigma^2)$ and exponentiated to terminal spot prices. Volatility held at the current VIX (flat).

**2. Two Risk Drivers — Bivariate Gaussian Model**
Joint estimation of (SPX log-return, VIX log-return) by MLE: vector mean and covariance matrix (with n-divisor MLE adjustment). Scenarios drawn from $\mathcal N_2(T\hat\mu, T\hat\Sigma)$ via `mvtnorm::rmvnorm`. Captures the empirical negative correlation between SPX and VIX returns — the natural delta/vega hedging effect: long-call portfolios are vega-positive, so simulated VIX spikes during equity drawdowns partially offset delta losses, producing a less-negative VaR/ES than the univariate model.

**3. Copula-Marginal Model — Student-t Marginals + Gaussian Copula**
Decoupled marginal and dependence structure:
- Fit Student-t marginals separately by MLE (`MASS::fitdistr`): df = 10 for SPX, df = 5 for VIX (heavier tails on VIX). Estimated location, scale, and converted scale to standard deviation via $\sigma = s\sqrt{df/(df-2)}$.
- Probability integral transform each series to uniforms via the Student-t CDF.
- Fit a Gaussian copula by MLE (`copula::fitCopula`) on the joint uniforms; recover correlation parameter.
- Simulate 5-day paths: draw daily uniform pairs from the fitted copula, transform back through Student-t inverse CDFs to log-returns, sum across 5 days for cumulative log-returns, exponentiate to terminal spot/VIX.

**4. Parametric Volatility Surface Model with VIX Shift**
Calibrated a parametric IV surface:

$$\sigma(m, \tau) = \alpha_1 + \alpha_2(m-1)^2 + \alpha_3(m-1)^3 + \alpha_4\sqrt{\tau}$$

where m = K/S is moneyness and τ is time-to-maturity in years. Fit the αs by minimizing absolute error to market IVs of OTM calls and OTM puts (avoiding ITM noise).

For pricing, the surface is dynamically anchored to the prevailing VIX:

$$\tilde\sigma(m,\tau) = \sigma(m,\tau) - (\alpha_1+\alpha_4) + \mathrm{VIX}_T$$

so that the ATM/short-tenor component tracks the VIX state. At T+5, $\mathrm{VIX}_{T+5}$ is taken from the bivariate Gaussian joint simulation, giving a coherent spot-VIX-surface scenario set.

**P&L Mechanics.** For each model the inception portfolio value is fixed at $V_0 = \sum_i \mathrm{BS}(S_0, K_i, T_i, \sigma_i, r_i)$. Each scenario $(S_{T+5}, \sigma_{T+5})$ revaluates the portfolio with the forward-curve discount factors at T+5 and the appropriate volatility input (flat VIX, simulated VIX, or VIX-shifted surface IV per option). 95% VaR and ES are computed from the empirical scenario P&L distribution (10,000 scenarios).

---

## Key Findings

- **Bivariate vs univariate**: Less-negative VaR/ES, since negative SPX–VIX correlation produces compensating vega gains during equity drawdowns.
- **Copula-marginal vs bivariate Gaussian**: Heavier-tailed Student-t marginals widen the loss tail; the Gaussian copula preserves linear dependence while allowing tail flexibility from the marginals.
- **Volatility surface vs flat VIX**: For this OTM-call portfolio, the calibrated surface assigns ~5–6 vol points above the ATM VIX (right side of the smile/skew). All four options have $\tilde\sigma - \mathrm{VIX}_T > 0$, so the surface model produces higher inception values and a wider, more negative P&L distribution than the bivariate Gaussian — capturing the smile premium that the flat-VIX model ignores.

---

## Key Files

| File | Description |
|------|-------------|
| [`main.Rmd`](main.Rmd) | Main R Markdown notebook — full analysis with all four models, plots, and commentary |
| [`main.html`](main.html) | Rendered HTML report with all outputs (recommended for review) |
| [`Functions/f_black_scholes.R`](Functions/f_black_scholes.R) | Black-Scholes pricing, forward rate construction from the term structure |
| [`Functions/f_pnl_plot.R`](Functions/f_pnl_plot.R) | P&L distribution computation (univariate, bivariate, copula, surface) and VaR/ES plotting |
| [`Functions/f_vol_surface_workflow.R`](Functions/f_vol_surface_workflow.R) | OTM data preparation, parametric surface fit, VIX-shift adjustment, surface-based P&L |
| [`Functions/f_copula_diagnostics_plot.R`](Functions/f_copula_diagnostics_plot.R) | Copula fit diagnostics: joint distribution, uniform pseudo-data, correlation checks |
| [`Outputs/`](Outputs/) | Generated plots: P&L distributions for all four models, fitted IV surface (3D HTML + 2D slices), copula diagnostics |
| [`docs/Instructions/`](docs/Instructions/) | Assignment description (MATH-60633A TP-02) |

**Data:** `Data/Raw/Market.rda` — pre-loaded `.rda` containing SPX series, VIX series, risk-free term structure, and the option chain (calls & puts).

---

## Tools & Libraries

R · `mvtnorm` (multivariate Gaussian sampling) · `copula` (Gaussian copula MLE, simulation) · `MASS` (Student-t MLE via `fitdistr`) · `xts` / `zoo` (time series) · `dplyr` / `purrr` (data manipulation) · `ggplot2` / `patchwork` (plots) · `plotly` (3D IV surface) · `here` (project paths) · `renv` (package reproducibility) · `knitr` (report generation)

---

## Setup

Open `TP2_Stats_Fin_Data.Rproj` in RStudio to activate `renv`, then restore packages:

```r
renv::restore()
```

Knit `main.Rmd` to regenerate the full HTML report (figures and tables under `Outputs/`).
