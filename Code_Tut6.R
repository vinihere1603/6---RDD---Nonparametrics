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



############################################################
# Tutorial 6 – Non-parametric Regression (R version)
# Data: Tutorial6_data2.dta
############################################################

library(haven)
library(dplyr)
library(np)          # for kernel regression and bandwidth selection
library(KernSmooth)  # for local polynomial (lpoly analogue)

#-----------------------------------------------------------
# Load data
#-----------------------------------------------------------

data_2 <- read_dta("C:/Users/kuehn/Meine Ablage/Uni/aktuelle Kurse/Impact Evaluation & Causal Inference/Tutorials/6 - RDD - Nonparametrics/Tutorial6_data2.dta")

# Ensure numeric where needed
data_2 <- data_2 %>%
  mutate(
    prestige  = as.numeric(prestige),
    income    = as.numeric(income),
    education = as.numeric(education),
    women     = as.numeric(women),
    type      = as.factor(type)
  )

prestige  <- data_2$prestige
income    <- data_2$income
education <- data_2$education
women     <- data_2$women
type_fac  <- data_2$type

############################################################
# A) Kernel Density Estimation (prestige)
############################################################

# A1: Histogram + kernel density (default Gaussian, default bw)
par(mfrow = c(1, 1))
hist(prestige,
     breaks = 20,
     prob   = TRUE,
     main   = "Prestige: Histogram + Kernel Density",
     xlab   = "Prestige")
lines(density(prestige), lwd = 2)

# A2: Kernel density with bandwidth = 1 (Gaussian)
dens_bw1 <- density(prestige, bw = 1, kernel = "gaussian")
plot(dens_bw1,
     main = "Kernel density of prestige (Gaussian, bw = 1)",
     xlab = "Prestige")

# A3: Kernel density with bandwidth = 20 (Gaussian)
dens_bw20 <- density(prestige, bw = 20, kernel = "gaussian")
plot(dens_bw20,
     main = "Kernel density of prestige (Gaussian, bw = 20)",
     xlab = "Prestige")

# A4: Kernel density with Gaussian kernel (let R choose bw)
dens_gauss <- density(prestige, kernel = "gaussian")
plot(dens_gauss,
     main = "Kernel density of prestige (Gaussian, automatic bw)",
     xlab = "Prestige")

# A5: Kernel density with rectangular (uniform) kernel
dens_rect <- density(prestige, kernel = "rectangular")
plot(dens_rect,
     main = "Kernel density of prestige (Rectangular, automatic bw)",
     xlab = "Prestige")

############################################################
# B) Kernel Regressions
############################################################

# Scatter + OLS line (unchanged)
plot(income, prestige,
     pch  = 16, cex = 0.6,
     xlab = "Income",
     ylab = "Prestige",
     main = "Scatter: prestige vs income")

ols_lin <- lm(prestige ~ income, data = data_2)
abline(ols_lin, col = "red", lwd = 2)

# Grid over income
grid_inc <- seq(min(income, na.rm = TRUE),
                max(income, na.rm = TRUE),
                length.out = 200)

## B3: local linear npreg with fixed bandwidth = 2000
mod_2000 <- npreg(
  xdat    = income,
  ydat    = prestige,
  regtype = "ll",
  bws     = 2000,                 # FIXED bandwidth
  exdat   = grid_inc              # evaluate on grid
)

pred_2000 <- mod_2000$mean

plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income", ylab = "Prestige",
     main = "Nonparametric regression: bw = 2000")
lines(grid_inc, pred_2000, col = "blue", lwd = 2)

## B5: fixed bandwidth = 5000
mod_5000 <- npreg(
  xdat    = income,
  ydat    = prestige,
  regtype = "ll",
  bws     = 5000,                 # FIXED bandwidth
  exdat   = grid_inc
)

pred_5000 <- mod_5000$mean

plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income", ylab = "Prestige",
     main = "Nonparametric regression: bw = 5000")
lines(grid_inc, pred_5000, col = "blue", lwd = 2)

## B6: automatic bandwidth via CV-LS (THIS is where npregbw belongs)
bw_auto_cv <- npregbw(
  xdat     = income,
  ydat     = prestige,
  regtype  = "ll",
  bwmethod = "cv.ls"
)

mod_auto_cv <- npreg(
  bws   = bw_auto_cv,
  exdat = grid_inc
)

pred_auto_cv <- mod_auto_cv$mean

plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income", ylab = "Prestige",
     main = "Nonparametric regression: automatic bw (CV-LS)")
lines(grid_inc, pred_auto_cv, col = "blue", lwd = 2)

## B7: automatic bandwidth via AIC-type criterion
bw_auto_aic <- npregbw(
  xdat     = income,
  ydat     = prestige,
  regtype  = "ll",
  bwmethod = "cv.aic"
)

mod_auto_aic <- npreg(
  bws   = bw_auto_aic,
  exdat = grid_inc
)

pred_auto_aic <- mod_auto_aic$mean

plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income", ylab = "Prestige",
     main = "Nonparametric regression: automatic bw (AIC-like)")
lines(grid_inc, pred_auto_aic, col = "blue", lwd = 2)

## B8: compare all four fits
par(mfrow = c(2, 2))

# (a) bw = 2000
plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income", ylab = "Prestige",
     main = "bw = 2000")
lines(grid_inc, pred_2000, col = "blue", lwd = 2)

# (b) bw = 5000
plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income", ylab = "Prestige",
     main = "bw = 5000")
lines(grid_inc, pred_5000, col = "blue", lwd = 2)

# (c) auto CV-LS
plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income", ylab = "Prestige",
     main = "Auto bw (CV-LS)")
lines(grid_inc, pred_auto_cv, col = "blue", lwd = 2)

# (d) auto AIC
plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income", ylab = "Prestige",
     main = "Auto bw (AIC-like)")
lines(grid_inc, pred_auto_aic, col = "blue", lwd = 2)

par(mfrow = c(1, 1))

############################################################
# B9: effect of income at medians of women & education
############################################################

# bandwidth selection (as before)
bw_multi <- npregbw(
  xdat    = data.frame(income    = income,
                       women     = women,
                       education = education),
  ydat    = prestige,
  regtype = "ll",
  bwmethod = "cv.ls"
)

# keep this if you want summary() etc.
mod_np_multi <- npreg(bws = bw_multi)
summary(mod_np_multi)

med_edu   <- median(education, na.rm = TRUE)
med_women <- median(women,     na.rm = TRUE)
med_inc   <- median(income,    na.rm = TRUE)

grid_income_multi <- data.frame(
  income    = seq(min(income, na.rm = TRUE),
                  max(income, na.rm = TRUE),
                  length.out = 200),
  women     = med_women,
  education = med_edu
)

# evaluate np-regression on this grid
mod_multi_B9 <- npreg(
  bws   = bw_multi,
  exdat = grid_income_multi
)

pred_multi <- mod_multi_B9$mean   # length == nrow(grid_income_multi)

plot(
  grid_income_multi$income,
  pred_multi,
  type = "l",
  xlab = "Income",
  ylab = "Predicted prestige",
  main = "NP regression: prestige ~ income + women + education\n(women & education at medians)"
)


############################################################
# B10: marginal effect of women (finite differences)
############################################################

grid_women <- seq(min(women, na.rm = TRUE),
                  max(women, na.rm = TRUE),
                  length.out = 50)

new_w_grid <- data.frame(
  income    = med_inc,
  women     = grid_women,
  education = med_edu
)

mod_multi_B10 <- npreg(
  bws   = bw_multi,
  exdat = new_w_grid
)

pred_w <- mod_multi_B10$mean   # length == length(grid_women)

# finite-difference approximation to dE[prestige]/dwomen
me_women <- c(diff(pred_w) / diff(grid_women), NA)

plot(
  grid_women[-length(grid_women)],
  me_women[-length(me_women)],
  type = "l",
  xlab = "Women (share in occupation)",
  ylab = "Approx. marginal effect on prestige",
  main = "Approx. marginal effect of women on prestige"
)


############################################################
# B11: prestige on income and type (factor)
############################################################

xdat_type <- data.frame(
  income = income,
  type   = type_fac
)

bw_type <- npregbw(
  xdat     = xdat_type,
  ydat     = prestige,
  regtype  = "ll",
  bwmethod = "cv.ls"
)

mod_type <- npreg(bws = bw_type)  # keep this if you want summary()
summary(mod_type)

############################################################
# B12: "marginal effects" of type at mean income
############################################################

mean_inc    <- mean(income, na.rm = TRUE)
type_levels <- levels(type_fac)

grid_type <- data.frame(
  income = rep(mean_inc, length(type_levels)),
  type   = factor(type_levels, levels = type_levels)
)

# Evaluate on this grid via exdat, not predict()
mod_type_B12 <- npreg(
  bws   = bw_type,
  exdat = grid_type
)

pred_type <- mod_type_B12$mean          # numeric, length == length(type_levels)

base_pred <- pred_type[1]
me_type   <- pred_type - base_pred      # also length == length(type_levels)

barplot(
  height    = me_type,
  names.arg = type_levels,
  xlab      = "Occupation type",
  ylab      = "Difference in predicted prestige vs baseline type",
  main      = "Type effects on prestige at mean income"
)

############################################################
# B13: Expected prestige vs income, by type
############################################################

grid_inc_type <- seq(
  min(income, na.rm = TRUE),
  max(income, na.rm = TRUE),
  length.out = 100
)

new_grid <- expand.grid(
  income = grid_inc_type,
  type   = type_levels
)

mod_type_B13 <- npreg(
  bws   = bw_type,
  exdat = new_grid
)

pred_inc_type <- mod_type_B13$mean   # length == nrow(new_grid)

par(mfrow = c(1, 1))
plot(NULL,
     xlim = range(grid_inc_type),
     ylim = range(pred_inc_type, na.rm = TRUE),
     xlab = "Income",
     ylab = "Predicted prestige",
     main = "Predicted prestige vs income, by occupation type")

cols <- seq_along(type_levels)

for (j in seq_along(type_levels)) {
  idx <- new_grid$type == type_levels[j]
  lines(grid_inc_type, pred_inc_type[idx], col = cols[j], lwd = 2)
}

legend("topleft",
       legend = type_levels,
       col    = cols,
       lwd    = 2,
       bty    = "n")


############################################################
# C) Local polynomial regression (lpoly analogue)
############################################################

# C14: local constant (degree 0) using locpoly, approximate "local average"
bw_loc <- dpill(income, prestige)    # data-driven bandwidth
fit_loc0 <- locpoly(income, prestige,
                    degree    = 0,
                    bandwidth = bw_loc,
                    gridsize  = 200)

plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income",
     ylab = "Prestige",
     main = "Local constant regression (degree 0)")
lines(fit_loc0$x, fit_loc0$y, col = "blue", lwd = 2)

# C15: local polynomial degree 3 (local cubic)
fit_loc3 <- locpoly(income, prestige,
                    degree    = 3,
                    bandwidth = bw_loc,
                    gridsize  = 200)

plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income",
     ylab = "Prestige",
     main = "Local polynomial regression (degree 3, bw = dpill)")
lines(fit_loc3$x, fit_loc3$y, col = "blue", lwd = 2)

# C16: local polynomial degree 3, bandwidth = 2000
fit_loc3_2000 <- locpoly(income, prestige,
                         degree    = 3,
                         bandwidth = 2000,
                         gridsize  = 200)

plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income",
     ylab = "Prestige",
     main = "Local polynomial regression (degree 3, bw = 2000)")
lines(fit_loc3_2000$x, fit_loc3_2000$y, col = "blue", lwd = 2)

# C17: Put local polynomial graphs plus some NP graphs in same frame
par(mfrow = c(2, 3))

# Graph from B3 (bw = 2000, npreg)
plot(income, prestige,
     pch = 16, cex = 0.5,
     main = "np: bw = 2000",
     xlab = "Income", ylab = "Prestige")
lines(grid_inc, pred_2000, col = "blue", lwd = 2)

# Graph from B5 (bw = 5000)
plot(income, prestige,
     pch = 16, cex = 0.5,
     main = "np: bw = 5000",
     xlab = "Income", ylab = "Prestige")
lines(grid_inc, pred_5000, col = "blue", lwd = 2)

# Local constant
plot(income, prestige,
     pch = 16, cex = 0.5,
     main = "locpoly: deg 0, bw = dpill",
     xlab = "Income", ylab = "Prestige")
lines(fit_loc0$x, fit_loc0$y, col = "blue", lwd = 2)

# Local cubic (dpill)
plot(income, prestige,
     pch = 16, cex = 0.5,
     main = "locpoly: deg 3, bw = dpill",
     xlab = "Income", ylab = "Prestige")
lines(fit_loc3$x, fit_loc3$y, col = "blue", lwd = 2)

# Local cubic (bw = 2000)
plot(income, prestige,
     pch = 16, cex = 0.5,
     main = "locpoly: deg 3, bw = 2000",
     xlab = "Income", ylab = "Prestige")
lines(fit_loc3_2000$x, fit_loc3_2000$y, col = "blue", lwd = 2)

par(mfrow = c(1, 1))

############################################################
# D) Locally weighted regression (lowess analogue)
############################################################

# D18: locally weighted regression of prestige on income
lowess_fit <- lowess(income, prestige, f = 2/3)

plot(income, prestige,
     pch = 16, cex = 0.6,
     xlab = "Income",
     ylab = "Prestige",
     main = "Lowess: locally weighted regression")
lines(lowess_fit, col = "blue", lwd = 2)