library(cmdstanr) # we use STAN
library(coda)
library(posterior)
library(ggplot2)
library(loo)

set.seed(42)
data <- read.csv("ex3/logistic_data.csv")
summary(data)
data <- list(
  N = nrow(data),
  x = data$x,
  y = data$y
)
sprintf("Number of observations: %d", data$N)

df <- as.data.frame(data)
p2 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.5, size = 2) +
  labs(title = "Data Plot",
       x = "Covariate (x)", y = "Outcome (y)") +
  theme_minimal()
print(p2)

# 1. Posterior formulation and implementation
mod1 <- cmdstan_model("ex3/logistic_regressionM1.stan")


# 2. MCMC fitting and diagnostics
fit1 <- mod1$sample(data=data,
    chains=4,
    parallel_chains=4,
    iter_warmup=1000,
    iter_sampling=10000,
    refresh = 0, # Silences the iteration progress updates
)
fit1$summary(c("beta0", "beta1"))
fit1$cmdstan_diagnose()

mcmc_full <- as_mcmc.list(fit1)
mcmc_chain <- mcmc_full[, c("beta0", "beta1")]
summary(mcmc_chain)

# trace plots
mcmc_mat <- as.matrix(mcmc_chain)
par(mfrow=c(1,2))
for (j in 1:ncol(mcmc_mat)) {
    plot(mcmc_mat[,j], type="l", main=paste("Trace theta", j), ylab="", xlab="")
}

# ESS
ess <- effectiveSize(mcmc_mat)
print(ess)

# autocorrelation
par(mfrow=c(1, 2)) 
for (j in 1:ncol(mcmc_mat)) {
    acf(mcmc_mat[,j], main=paste("ACF for", colnames(mcmc_mat)[j]))
}


# 3. Posterior inference
fit1$summary(c("beta0", "beta1"))

# plot of the marginal posterior distributions
par(mfrow=c(1,2))
for (j in 1:ncol(mcmc_mat)) {
    hist(mcmc_mat[,j], breaks=30, main=paste("Posterior of", colnames(mcmc_mat)[j]), xlab="", freq=FALSE)
}

# plots of the joint posterior distribution
par(mfrow=c(1,1))
plot(mcmc_mat[,1], mcmc_mat[,2], main="Joint Posterior of beta0 and beta1", xlab="beta0", ylab="beta1", pch=16, cex=0.5)
par(mfrow=c(1,1)) 
df <- as_draws_df(fit1$draws())
vars <- c("beta0", "beta1")
cor(as.matrix(df[, vars]))
ggplot(df, aes(x = beta0, y = beta1)) + geom_density_2d_filled() + theme_minimal()


# 4. Posterior predictive inference
x_new <- 1.0
beta0_draws <- mcmc_mat[, "beta0"]
beta1_draws <- mcmc_mat[, "beta1"]

# plogis() is the inverse logit function
pi_new <- plogis(beta0_draws + beta1_draws * x_new)

# posterior mean and 95% credible interval for pi_new
pi_mean <- mean(pi_new)
pi_ci <- quantile(pi_new, probs = c(0.025, 0.975))

sprintf("Posterior mean of pi_new: %.4f", pi_mean)
sprintf("95%% Credible Interval for pi_new: [%.4f, %.4f]", pi_ci[1], pi_ci[2])

par(mfrow=c(1,1))
hist(pi_new, breaks=50, col="steelblue", border="white", 
     main=expression(paste("Posterior Distribution of ", pi[new])), 
     xlab=expression(pi[new]), freq=FALSE)
density_pi_new <- density(pi_new)
lines(density_pi_new, col="darkred", lwd=2)


# 5. Model comparison
mod0 <- cmdstan_model("ex3/logistic_regressionM0.stan")

# we fit and sample from M0
fit0 <- mod0$sample(data=data,
    chains=4,
    parallel_chains=4,
    iter_warmup=1000,
    iter_sampling=10000,
    refresh = 0, 
)
fit0$summary(c("beta0"))
fit0$cmdstan_diagnose()

# comparison using LOO
log_lik_M1 <- fit1$draws("log_lik")
log_lik_M0 <- fit0$draws("log_lik")

# compute LOO for both models
loo_M1 <- loo(log_lik_M1)
loo_M0 <- loo(log_lik_M0)
comparison <- loo_compare(loo_M0, loo_M1)
print(comparison)