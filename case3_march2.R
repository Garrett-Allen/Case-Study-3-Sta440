###Start with a clean space
rm(list = ls())

###Load packages
library(womblR)
library(rstan)
library(loo)
library(LaplacesDemon) # to sample from halt-t

###Simulate data 
set.seed(2)
T <- 15
n <- 54
N <- n * T
W <- womblR::HFAII_QueenHF
blind_spot <- c(26, 35) # define blind spot
# W <- HFAII_QueenHF[-blind_spot, -blind_spot] # visual field adjacency matrix
W <- HFAII_QueenHF # visual field adjacency matrix
D <- diag(apply(W, 1, sum))
# Q <- D - 0.99 * W
Q <- 0.99 * (D - W) + 0.01 * diag(n)
precision0 <- 1.4 * Q
sigma0 <- solve(precision0)
beta0 <- 26 + chol(sigma0) %*% matrix(rnorm(n), ncol = 1)
precision1 <- 0.5 * Q
sigma1 <- solve(precision1)
beta1 <- -6 + chol(sigma1) %*% matrix(rnorm(n), ncol = 1)
sigma <- rhalft(n, scale = 1.5, nu = 3)
y <- matrix(nrow = n, ncol = T)
time <- seq(0, 1, length.out = T)
for (i in 1:n) {
  for (t in 1:T) {
    mu <- beta0[i] + beta1[i] * time[t]
    y[i, t] <- rnorm(1, mu, sigma[i])
  }
}  

###Create training data
T_train <- 9
dat <- data.frame(y = as.numeric(y[, 1:T_train]), location = rep(1:n, T_train), time = rep(time[1:T_train], each = n))

###Plot the data
pdf("example.pdf", height = 5, width = 5)
PlotVfTimeSeries(Y = dat$y,
                 Location = dat$location,
                 Time = dat$time * 365, line.reg = FALSE)
dev.off()

###Define stan data list
blind_spot <- c(26, 35)
dat <- dat[!dat$location %in% blind_spot, ]
W <- W[-blind_spot, -blind_spot]
n_train <- 52
dat$location <- as.numeric(as.factor(dat$location)) # reorder the locations (1-52)
stan_data <- list(
  y = dat$y,
  x = dat$time,
  s = dat$location,
  N = T_train * n_train,
  n = n_train,
  T = T_train,
  W = W
)

###Compile spatial model
model_compiled <- stan_model("spatial.stan")
fit <- sampling(model_compiled, data = stan_data, chains = 5, iter = 5000, cores = 5)
saveRDS(fit, file = "fit.rds")

###Compile no spatial model
model_compiled_plr <- stan_model("plr.stan")
fit_plr <- sampling(model_compiled_plr, data = stan_data, chains = 5, iter = 5000, cores = 5)
saveRDS(fit_plr, file = "fit_plr.rds")

###Check summaries
print(fit, pars = c("beta0", "beta1", "sigma", "tau0", "rho0", "tau1", "rho1"), digits_summary = 5)

###Check traceplots
pdf("traceplots.pdf", height = 5, width = 9)
rstan::traceplot(fit, pars = c("beta0", "beta1", "tau0", "rho0", "tau1", "rho1"))
dev.off()

###Model fit
log_lik <- loo::extract_log_lik(fit)
waic <- loo::waic(log_lik)
waic$estimates[3, 1]
log_lik_plr <- loo::extract_log_lik(fit_plr)
waic_plr <- loo::waic(log_lik_plr)
waic_plr$estimates[3, 1]

###Get posterior predictive distribution for spatial version, and compute residuals
params <- extract(fit)
beta0 <- params$beta0
beta1 <- params$beta1
beta0_vec <- params$beta0_vec
beta1_vec <- params$beta1_vec
sigma <- params$sigma
n_sims <- length(beta0)
pred <- array(dim = c(n_sims, n_train, T))
for (sim in 1:n_sims) {
  for (i in 1:n_train) {
    for (t in 1:T) {
      mu <- (beta0[sim] + beta0_vec[sim, i]) + (beta1[sim] + beta1_vec[sim, i]) * time[t]
      pred[sim, i, t] <- rnorm(1, mu, sigma[sim, i])
    }
  }
}
ppd_mean <- apply(pred, c(2, 3), mean)
ppd_sd <- apply(pred, c(2, 3), sd)
y_train <- y[-c(26, 35), ]
resids <- y_train - ppd_mean

###Posterior predictive distribution for PLR model, and compute residuals
params <- extract(fit_plr)
beta0 <- params$beta0
beta1 <- params$beta1
sigma <- params$sigma
n_sims <- nrow(beta0)
pred <- array(dim = c(n_sims, n_train, T))
for (sim in 1:n_sims) {
  for (i in 1:n_train) {
    for (t in 1:T) {
      mu <- beta0[sim, i] + beta1[sim, i] * time[t]
      pred[sim, i, t] <- rnorm(1, mu, sigma[sim, i])
    }
  }
}  
ppd_mean_plr <- apply(pred, c(2, 3), mean)
ppd_sd_plr <- apply(pred, c(2, 3), sd)
resids_plr <- y_train - ppd_mean_plr

###MSE plot
pdf("prediction.pdf", height = 5, width = 7)
plot(time, apply(abs(resids)^2, 2, mean), ylim = c(0, 15), xlab = "Time", ylab = "MSE", pch = 15, col = "black")
points(time, apply(abs(resids_plr)^2, 2, mean), pch = 15, col = "red")
legend("topleft", legend = c("Spatial model", "PLR"), pch = 15, col = c("black", "red"), bty = "n")
abline(v = mean(time[9:10])) # cutoff for future predictions
dev.off()


###Look at distribution of slopes and intercepts
params <- extract(fit)
beta0 <- mean(params$beta0)
beta1 <- mean(params$beta1)
beta0_vec <- apply(params$beta0_vec, 2, mean)
beta1_vec <- apply(params$beta1_vec, 2, mean)
sigma <- apply(params$sigma, 2, mean)
intercepts <- beta0 + beta0_vec
slopes <- beta1 + beta1_vec
params <- extract(fit_plr)
intercepts_plr <- apply(params$beta0, 2, mean)
slopes_plr <- apply(params$beta1, 2, mean)
pdf("slopes.pdf", height = 6, width = 6)
PlotSensitivity(Y = slopes,
                main = "Spatial Slopes",
                legend.lab = "DLS (dB) / year", legend.round = 2,
                bins = 250, border = FALSE, zlim = c(-20, 10))
dev.off()
pdf("slopes_plr.pdf", height = 6, width = 6)
PlotSensitivity(Y = slopes_plr,
                main = "PLR Slopes",
                legend.lab = "DLS (dB)", legend.round = 2,
                bins = 250, border = FALSE, zlim = c(-20, 10))
dev.off()
pdf("intercepts.pdf", height = 6, width = 6)
PlotSensitivity(Y = intercepts,
                main = "Spatial Intercepts",
                legend.lab = "DLS (dB)", legend.round = 2,
                bins = 250, border = FALSE, zlim = c(20, 34))
dev.off()
pdf("intercepts_plr.pdf", height = 6, width = 6)
PlotSensitivity(Y = intercepts_plr,
                main = "PLR Intercepts",
                legend.lab = "DLS (dB)", legend.round = 2,
                bins = 250, border = FALSE, zlim = c(20, 34))
dev.off()
