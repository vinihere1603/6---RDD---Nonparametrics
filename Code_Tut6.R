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



### 2. Nonparametric regression  ------------------------------------------------

library(haven)
library(dplyr)
library(KernSmooth)
library(np)

## Load data --------------------------------------------------------------------
data_2 <- read_dta(
  "C:/Users/kuehn/Meine Ablage/Uni/aktuelle Kurse/Impact Evaluation & Causal Inference/Tutorials/6 - RDD - Nonparametrics/Tutorial6_data2.dta"
)

data_2 <- data_2 %>%
  mutate(type = factor(type))

# Shorthands
prestige  <- data_2$prestige
income    <- data_2$income
women     <- data_2$women
education <- data_2$education
type      <- data_2$type


## (1) Kernel density of prestige ----------------------------------------------

par(mfrow = c(1, 1))

dens_prestige <- density(prestige, na.rm = TRUE, kernel = "gaussian")

plot(
  dens_prestige,
  main = "Kernel density of prestige",
  xlab = "Prestige"
)


## (2) Nonparametric regression prestige ~ income -------------------------------

# Scatter + local polynomial + LOWESS
plot(
  income, prestige,
  pch = 16, cex = 0.6,
  xlab = "Income",
  ylab = "Prestige",
  main = "Prestige vs Income"
)

# Local polynomial (KernSmooth, degree 1)
bw_lp <- sd(income, na.rm = TRUE) / 4
lp_inc <- locpoly(income, prestige, degree = 1, bandwidth = bw_lp)
lines(lp_inc, lwd = 2)

# LOWESS
lw_inc <- lowess(income, prestige, f = 0.5)
lines(lw_inc, lwd = 2, lty = 2)

legend(
  "topleft",
  legend = c("local poly", "LOWESS"),
  lwd = 2,
  lty = c(1, 2),
  bty = "n"
)


## (3) npreg: prestige ~ income (univariate) ------------------------------------

bw_inc_np <- npregbw(
  prestige ~ income,
  data     = data_2,
  regtype  = "lc",
  bwmethod = "cv.ls"
)
bw_inc_np

np_inc <- npreg(bws = bw_inc_np)
summary(np_inc)

# Plot fitted CEF
inc_grid <- seq(
  from = min(income, na.rm = TRUE),
  to   = max(income, na.rm = TRUE),
  length.out = 200
)
pred_inc <- predict(np_inc, newdata = data.frame(income = inc_grid))

plot(
  income, prestige,
  pch = 16, cex = 0.5, col = "grey70",
  xlab = "Income",
  ylab = "Prestige",
  main = "npreg: E[Prestige | Income]"
)
lines(inc_grid, pred_inc, lwd = 2)


## (4) Multivariate npreg: prestige ~ income + women + education ---------------

bw_multi <- npregbw(
  prestige ~ income + women + education,
  data     = data_2,
  regtype  = "lc",
  bwmethod = "cv.ls"
)
bw_multi

np_multi <- npreg(bws = bw_multi)
summary(np_multi)

# Marginal effect of women: vary women, fix income and education at means
inc_bar <- mean(income, na.rm = TRUE)
edu_bar <- mean(education, na.rm = TRUE)

grid_women <- seq(
  from = min(women, na.rm = TRUE),
  to   = max(women, na.rm = TRUE),
  length.out = 50
)

new_multi <- data.frame(
  income    = rep(inc_bar, length(grid_women)),
  women     = grid_women,
  education = rep(edu_bar, length(grid_women))
)

pred_w <- predict(np_multi, newdata = new_multi)

plot(
  grid_women, pred_w,
  type = "l",
  xlab = "Women share in occupation",
  ylab = "E[Prestige | income, women, education]",
  main = "Effect of women share (income, education at means)"
)


## (5) npreg with discrete covariate: prestige ~ income + type ------------------

bw_type <- npregbw(
  prestige ~ income + type,
  data     = data_2,
  regtype  = "lc",
  bwmethod = "cv.ls"
)
bw_type

np_type <- npreg(bws = bw_type)
summary(np_type)

# Predicted prestige at median income by occupation type
type_levels <- levels(type)
inc_med     <- median(income, na.rm = TRUE)

pred_type <- sapply(type_levels, function(tk) {
  newdata <- data.frame(
    income = inc_med,
    type   = factor(tk, levels = type_levels)
  )
  predict(np_type, newdata = newdata)
})

pred_type_df <- data.frame(
  type         = type_levels,
  prestige_hat = pred_type
)
print(pred_type_df)


## (6) CEF of prestige vs income, separate curves by type ----------------------

income_grid <- seq(
  from = min(income, na.rm = TRUE),
  to   = max(income, na.rm = TRUE),
  length.out = 100
)

# Initialize empty plot
plot(
  NA,
  xlim = range(income_grid),
  ylim = range(prestige, na.rm = TRUE),
  xlab = "Income",
  ylab = "E[Prestige | Income, Type]",
  main = "npreg CEF: prestige ~ income + type"
)

cols <- seq_along(type_levels)

for (j in seq_along(type_levels)) {
  newdat_j <- data.frame(
    income = income_grid,
    type   = factor(type_levels[j], levels = type_levels)
  )
  yhat_j <- predict(np_type, newdata = newdat_j)
  lines(income_grid, yhat_j, lwd = 2, col = cols[j])
}

legend(
  "topleft",
  legend = type_levels,
  col    = cols,
  lwd    = 2,
  bty    = "n",
  title  = "Occupation type"
)
