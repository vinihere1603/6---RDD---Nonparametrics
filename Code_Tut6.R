### RDD and Nonparametrics

##1. RDD

#1.
# - fuzzy RDD: elegebility for longer unemeployment benefits is deterministic but uptake isn't
#   DiD between workers <50 and >50 before and after treatment

#2./3.
library(haven)
data_1  <- read_dta("C:/Users/kuehn/Meine Ablage/Uni/aktuelle Kurse/Impact Evaluation & Causal Inference/Tutorials/6 - RDD - Nonparametrics/Tutorial6_data1.dta")

library(np)
h    <- 2                 # bandwidth
a0   <- 50                # evaluation age (cutoff)
z    <- (data_1$age - a0) / h   # standardized distance

## Kernels -----------------------------------------------------------
# Triangular: K(z) = (1 - |z|) if |z|<1 else 0
K_tri <- function(z) ifelse(abs(z) < 1, 1 - abs(z), 0)

# Epanechnikov as in sheet: K(z) = (3/4)*(1 - z^2/5)/sqrt(5) if |z|<sqrt(5) else 0
K_epa <- function(z) {
  out  <- numeric(length(z))
  idx  <- abs(z) < sqrt(5)
  out[idx] <- (3/4) * (1 - z[idx]^2 / 5) / sqrt(5)
  out
}

# Normal kernel: standard normal density
K_norm <- function(z) dnorm(z)

## Un-normalized kernel values
k_tri  <- K_tri(z)
k_epa  <- K_epa(z)
k_norm <- K_norm(z)

## Normalized weights (sum to 1 for each kernel)
w_tri  <- k_tri  / sum(k_tri)
w_epa  <- k_epa  / sum(k_epa)
w_norm <- k_norm / sum(k_norm)

## Combine in data frame
ker_df <- data.frame(
  age    = data_1$age,
  w_tri  = w_tri,
  w_epa  = w_epa,
  w_norm = w_norm
)

## Plots -------------------------------------------------------------

# Base R: three panels
par(mfrow = c(1, 3))

plot(ker_df$age, ker_df$w_tri,
     main = "Triangular kernel",
     xlab = "Age", ylab = "Weight", pch = 16, cex = 0.6)

plot(ker_df$age, ker_df$w_epa,
     main = "Epanechnikov kernel",
     xlab = "Age", ylab = "Weight", pch = 16, cex = 0.6)

plot(ker_df$age, ker_df$w_norm,
     main = "Normal kernel",
     xlab = "Age", ylab = "Weight", pch = 16, cex = 0.6)

par(mfrow = c(1, 1))


#4.
library(dplyr)
data_1 <- data_1 %>%
  mutate(age_dev = age - a0)

#5.
# 1. Standardized distance and triangular kernel weights
age_dev <- data_1$age_dev          # deviation from cutoff
z       <- age_dev / h
k_tri   <- K_tri(z)                # unnormalized weights

# 2. Local linear regressions, right and left of cutoff

# Right side: age_dev >= 0  ->  y^+
fit_plus <- lm(
  unemployment_duration ~ age_dev,
  data    = data_1,
  weights = k_tri,
  subset  = age_dev >= 0
)
y_plus <- coef(fit_plus)[1]

# Left side: age_dev < 0  ->  y^-
fit_minus <- lm(
  unemployment_duration ~ age_dev,
  data    = data_1,
  weights = k_tri,
  subset  = age_dev < 0
)
y_minus <- coef(fit_minus)[1]

y_plus
y_minus

#6. CEF plot showing the discontinuity at the cutoff

cutoff <- a0

# Restrict to the bandwidth window: |age_dev| <= h
# This is the region where the triangular kernel actually has support.
grid_left_dev  <- seq(-h, 0, length.out = 100)
grid_right_dev <- seq(0,  h, length.out = 100)

# Predicted CEF from the two local-linear fits (piecewise around cutoff)
yhat_left  <- coef(fit_minus)[1] + coef(fit_minus)[2] * grid_left_dev
yhat_right <- coef(fit_plus)[1]  + coef(fit_plus)[2]  * grid_right_dev

# Back to age
age_left  <- cutoff + grid_left_dev
age_right <- cutoff + grid_right_dev

# Common y-limits for visibility of both sides
ylims <- range(c(yhat_left, yhat_right), na.rm = TRUE)

# Plot only the CEFs, no scatter
plot(
  age_left, yhat_left,
  type = "l",
  xlab = "Age",
  ylab = "E[Unemployment Duration | Age]",
  main = "RDD CEF: Local-Linear Triangular Kernel",
  xlim = cutoff + c(-h, h),
  ylim = ylims
)

lines(age_right, yhat_right)

abline(v = cutoff, lty = 2)  # mark cutoff


#7. Estimate the treatment effect using the RD estimator ----------------------

# Define treatment indicator as in sharp RDD: 1 if age >= cutoff, 0 otherwise
data_1 <- data_1 %>%
  mutate(treat = ifelse(age >= cutoff, 1, 0))

# RD regression with same local-linear spec and triangular weights:
fit_rd <- lm(
  unemployment_duration ~ treat + age_dev + treat:age_dev,
  data    = data_1,
  weights = k_tri
)

summary(fit_rd)

# RD estimator (treatment effect at cutoff)
tau_hat_rd <- coef(fit_rd)["treat"]
tau_hat_rd

# Check consistency with difference in limits from #5
tau_hat_limits <- y_plus - y_minus
tau_hat_limits


#9. Reproduce #5 and #6 using rdrobust (analogue of Stata's rd/rdrobust) -----

# install.packages("rdrobust") # run once
library(rdrobust)

# Local-linear, triangular kernel, manual bandwidth h on both sides
rd_out <- rdrobust(
  y       = data_1$unemployment_duration,
  x       = data_1$age,
  c       = cutoff,
  h       = c(h, h),              # left and right bandwidths
  p       = 1,                    # local linear
  kernel  = "triangular"
)

summary(rd_out)

# rdplot: CEF on both sides with the same kernel/bandwidth choice
rdplot(
  y      = data_1$unemployment_duration,
  x      = data_1$age,
  c      = cutoff,
  kernel = "triangular",
  p      = 1,
  h      = h,
  x.label = "Age",
  y.label = "Unemployment Duration",
  title   = "RDD CEF: rdrobust (local-linear, triangular)"
)


#10. Compute optimal bandwidth ----------------------------------------------

bw_opt <- rdbwselect(
  y       = data_1$unemployment_duration,
  x       = data_1$age,
  c       = cutoff,
  kernel  = "triangular",
  p       = 1
)

bw_opt
bw_opt$bws  # shows: row 'mserd', columns 'h (left)', 'h (right)', ...

# Extract MSE-optimal h on each side
h_opt_left  <- bw_opt$bws["mserd", "h (left)"]
h_opt_right <- bw_opt$bws["mserd", "h (right)"]

h_opt_left
h_opt_right

h_opt <- h_opt_left


#11. Estimate effects for different bandwidths and plot vs bandwidth ---------

bw_grid <- seq(0.5, 10, by = 0.5)

tau_bw <- sapply(bw_grid, function(bw) {
  rd_tmp <- rdrobust(
    y       = data_1$unemployment_duration,
    x       = data_1$age,
    c       = cutoff,
    h       = c(bw, bw),        # symmetric bandwidth L/R
    p       = 1,
    kernel  = "triangular"
  )
  rd_tmp$coef[1]   # treatment effect at cutoff
})

plot(
  bw_grid, tau_bw,
  type = "b",
  xlab = "Bandwidth",
  ylab = "RD estimate (tau_hat)",
  main = "RD estimates vs bandwidth (triangular, local-linear)"
)
abline(h = 0,       lty = 2)
abline(v = h,       lty = 3)  # your original manual bandwidth (h <- 2)
abline(v = h_opt,   lty = 3)  # MSE-optimal bandwidth from rdbwselect


#12. Placebo test with available variables ----------------------------------
# With only {age, unemployment_duration, treat}, you cannot do covariate
# balance tests or outcome placebo on other variables.
#
# But you CAN do a placebo cutoff test:
# - Choose false cutoffs where no policy change occurs, e.g. 45 or 55.
# - Re-estimate the RD effect at those fake cutoffs; you should find effects
#   close to zero if the design is credible and there is no underlying trend
#   "masquerading" as a treatment jump.

# Example: placebo cutoff at 45
cutoff_placebo <- 47

rd_placebo <- rdrobust(
  y       = data_1$unemployment_duration,
  x       = data_1$age,
  c       = cutoff_placebo,
  kernel  = "triangular",
  p       = 1
)

summary(rd_placebo)

# If the RD estimate at 47 is statistically small and imprecise, that supports
# the claim that the jump at 50 is due to the program rather than a smooth trend.


#13. Other tests with more variables (conceptual, no code) -------------------
# If you had additional variables, two natural tests:
#
# (1) Covariate balance tests at the cutoff:
#     - Use predetermined covariates (education, pre-program earnings, gender,
#       etc.) as outcomes in the same RD specification.
#     - Check that there is no significant jump in these covariates at 50.
#       Any discontinuity would signal selection or manipulation.
#
# (2) Density test for manipulation of the running variable:
#     - Use a McCrary-type density test (e.g. in rdrobust: rddensity) on age
#       around the cutoff.
#     - A significant jump in the density at 50 would indicate manipulation
#       of the running variable (people bunching just below/above the cutoff).
#
# Other possibilities with richer data:
# - Placebo outcomes that should not respond to the program (e.g. health or
#   non-labor outcomes unaffected by the treatment).
# - Heterogeneity analyses (different subgroups) to see if effects behave in
#   economically plausible ways.



# ============================================================
# Tutorial 6 – Part 1: Non-parametric Regression
# Data 2: prestige, income, education, women, census, type
# ============================================================

# Packages
library(haven)   # read_dta
library(np)      # nonparametrics

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
data_2 <- read_dta("C:/Users/kuehn/Meine Ablage/Uni/aktuelle Kurse/Impact Evaluation & Causal Inference/Tutorials/6 - RDD - Nonparametrics/Tutorial6_data2.dta")

# Quick check
str(data_2)
summary(data_2$prestige)
summary(data_2$income)

# ============================================================
# Part 1A – Kernel density estimation for prestige
# ============================================================

prestige <- data_2$prestige

# 1) Histogram + default kernel density estimate on same graph
par(mfrow = c(1, 1))
hist(
  prestige,
  breaks = 20,
  prob = TRUE,
  main = "Prestige: Histogram and default kernel density",
  xlab = "Prestige"
)
lines(
  density(prestige, na.rm = TRUE),
  lwd = 2
)

# 2) Kernel density with bandwidth = 1
dens_bw1 <- density(prestige, bw = 1, na.rm = TRUE)

plot(
  dens_bw1,
  main = "Prestige: Kernel density, bw = 1",
  xlab = "Prestige"
)

# 3) Kernel density with bandwidth = 20
dens_bw20 <- density(prestige, bw = 20, na.rm = TRUE)

plot(
  dens_bw20,
  main = "Prestige: Kernel density, bw = 20",
  xlab = "Prestige"
)

# 4) Kernel density with Gaussian kernel (explicit)
#    (base::density is Gaussian by default, this just makes it explicit)
dens_gauss <- density(prestige, kernel = "gaussian", na.rm = TRUE)

plot(
  dens_gauss,
  main = "Prestige: Kernel density, Gaussian kernel",
  xlab = "Prestige"
)

# 5) Kernel density with rectangular (uniform) kernel using np::npudens
#    np chooses bandwidth by cross-validation, we only change the kernel type.

# Gaussian version in np (for comparison)
dens_np_gauss <- npudens(
  ~ prestige,
  data     = data_2,
  ckertype = "gaussian"
)

plot(
  dens_np_gauss,
  main = "Prestige: np kernel density, Gaussian kernel",
  xlab = "Prestige"
)

# Uniform (rectangular) kernel in np
dens_np_uniform <- npudens(
  ~ prestige,
  data     = data_2,
  ckertype = "uniform"
)

plot(
  dens_np_uniform,
  main = "Prestige: np kernel density, uniform (rectangular) kernel",
  xlab = "Prestige"
)

# ============================================================
# Part 1B – Kernel regressions: prestige on income
# ============================================================

income   <- data_2$income
prestige <- data_2$prestige

# ------------------------------------------------------------
# B1) Parametric vs nonparametric regression
# ------------------------------------------------------------

# Parametric (OLS) polynomial in income if you want, here linear:
mod_lin <- lm(prestige ~ income, data = data_2)
summary(mod_lin)

# Nonparametric local linear regression (automatic bandwidth, cv.ls)
bw_auto <- npregbw(prestige ~ income, data = data_2, regtype = "ll")
bw_auto
mod_np_auto <- npreg(bws = bw_auto, data = data_2)

# Evaluation grid: 50 points
grid_50 <- data.frame(
  income = seq(
    min(income, na.rm = TRUE),
    max(income, na.rm = TRUE),
    length.out = 50
  )
)

pred_lin_50    <- predict(mod_lin,     newdata = grid_50)
pred_np_auto_50 <- predict(mod_np_auto, newdata = grid_50)

# Plot: scatter + OLS + nonparametric
par(mfrow = c(1, 1))
plot(
  income, prestige,
  pch  = 16,
  cex  = 0.6,
  xlab = "Income",
  ylab = "Prestige",
  main = "Prestige vs income: OLS vs nonparametric (50 grid points)"
)
lines(grid_50$income, pred_lin_50,      lwd = 2, lty = 2)
lines(grid_50$income, pred_np_auto_50,  lwd = 2)

# ------------------------------------------------------------
# B2) Same comparison with a finer grid: 200 evaluation points
# ------------------------------------------------------------

grid_200 <- data.frame(
  income = seq(
    min(income, na.rm = TRUE),
    max(income, na.rm = TRUE),
    length.out = 200
  )
)

pred_lin_200    <- predict(mod_lin,     newdata = grid_200)
pred_np_auto_200 <- predict(mod_np_auto, newdata = grid_200)

plot(
  income, prestige,
  pch  = 16,
  cex  = 0.6,
  xlab = "Income",
  ylab = "Prestige",
  main = "Prestige vs income: OLS vs nonparametric (200 grid points)"
)
lines(grid_200$income, pred_lin_200,      lwd = 2, lty = 2)
lines(grid_200$income, pred_np_auto_200,  lwd = 2)

# ------------------------------------------------------------
# B3) Nonparametric regression with fixed bandwidth = 2000,
#     bootstrapped standard errors (100 replications)
# ------------------------------------------------------------

mod_np_bw2000 <- npreg(
  prestige ~ income,
  data        = data_2,
  regtype     = "ll",
  bws         = 2000,          # fixed bandwidth
  boot.method = "wild",        # or "residual"
  boot.num    = 100
)

summary(mod_np_bw2000)

# ------------------------------------------------------------
# B4) Plot with bw = 2000
# ------------------------------------------------------------

grid_bw <- grid_200

pred_np_bw2000 <- predict(mod_np_bw2000, newdata = grid_bw)

plot(
  income, prestige,
  pch  = 16,
  cex  = 0.6,
  xlab = "Income",
  ylab = "Prestige",
  main = "Nonparametric regression: bw = 2000"
)
lines(grid_bw$income, pred_np_bw2000, lwd = 2)

# ------------------------------------------------------------
# B5) Nonparametric regression with bw = 5000
# ------------------------------------------------------------

mod_np_bw5000 <- npreg(
  prestige ~ income,
  data        = data_2,
  regtype     = "ll",
  bws         = 5000,
  boot.method = "wild",
  boot.num    = 100
)

summary(mod_np_bw5000)

pred_np_bw5000 <- predict(mod_np_bw5000, newdata = grid_bw)

plot(
  income, prestige,
  pch  = 16,
  cex  = 0.6,
  xlab = "Income",
  ylab = "Prestige",
  main = "Nonparametric regression: bw = 5000"
)
lines(grid_bw$income, pred_np_bw5000, lwd = 2)

# ------------------------------------------------------------
# B6) Automatic bandwidth (cv.ls) – already computed as bw_auto
#     This is the standard cross-validated least-squares bw.
# ------------------------------------------------------------

bw_auto
summary(mod_np_auto)

pred_np_auto_bw <- predict(mod_np_auto, newdata = grid_bw)

plot(
  income, prestige,
  pch  = 16,
  cex  = 0.6,
  xlab = "Income",
  ylab = "Prestige",
  main = "Nonparametric regression: automatic bw (cv.ls)"
)
lines(grid_bw$income, pred_np_auto_bw, lwd = 2)

# ------------------------------------------------------------
# B7) AIC-based bandwidth (analogue to Stata’s imaic)
# ------------------------------------------------------------

bw_aic <- npregbw(
  prestige ~ income,
  data     = data_2,
  regtype  = "ll",
  bwmethod = "cv.aic"
)

bw_aic
mod_np_aic <- npreg(bws = bw_aic, data = data_2)

pred_np_aic <- predict(mod_np_aic, newdata = grid_bw)

plot(
  income, prestige,
  pch  = 16,
  cex  = 0.6,
  xlab = "Income",
  ylab = "Prestige",
  main = "Nonparametric regression: AIC-based bw (cv.aic)"
)
lines(grid_bw$income, pred_np_aic, lwd = 2)

# ------------------------------------------------------------
# B8) Four-panel comparison of different bandwidth choices
# ------------------------------------------------------------

par(mfrow = c(2, 2))

# bw = 2000
plot(
  income, prestige,
  pch  = 16, cex = 0.5,
  xlab = "Income", ylab = "Prestige",
  main = "bw = 2000"
)
lines(grid_bw$income, pred_np_bw2000, lwd = 2)

# bw = 5000
plot(
  income, prestige,
  pch  = 16, cex = 0.5,
  xlab = "Income", ylab = "Prestige",
  main = "bw = 5000"
)
lines(grid_bw$income, pred_np_bw5000, lwd = 2)

# automatic (cv.ls)
plot(
  income, prestige,
  pch  = 16, cex = 0.5,
  xlab = "Income", ylab = "Prestige",
  main = "automatic bw (cv.ls)"
)
lines(grid_bw$income, pred_np_auto_bw, lwd = 2)

# AIC-based (cv.aic)
plot(
  income, prestige,
  pch  = 16, cex = 0.5,
  xlab = "Income", ylab = "Prestige",
  main = "AIC-based bw (cv.aic)"
)
lines(grid_bw$income, pred_np_aic, lwd = 2)

par(mfrow = c(1, 1))  # reset

# ------------------------------------------------------------
# B9) Multivariate nonparametric regression:
#     prestige on income, education, and women
# ------------------------------------------------------------

# Treat women as a factor (discrete regressor)
data_2$women_factor <- factor(data_2$women)

bw_multi <- npregbw(
  prestige ~ income + education + women_factor,
  data     = data_2,
  regtype  = "ll",
  bwmethod = "cv.aic"
)

bw_multi
mod_np_multi <- npreg(bws = bw_multi, data = data_2)
summary(mod_np_multi)

# Example: partial effect of income at median education and women = 0
med_edu   <- median(data_2$education, na.rm = TRUE)
base_w    <- levels(data_2$women_factor)[1]

grid_income_multi <- data.frame(
  income        = seq(
    min(income, na.rm = TRUE),
    max(income, na.rm = TRUE),
    length.out = 200
  ),
  education     = med_edu,
  women_factor  = base_w
)

pred_multi <- predict(mod_np_multi, newdata = grid_income_multi)

plot(
  grid_income_multi$income,
  pred_multi,
  type = "l",
  xlab = "Income",
  ylab = "Predicted prestige",
  main = "Nonparametric regression: prestige ~ income + education + women\n(education at median, women = base category)"
)

