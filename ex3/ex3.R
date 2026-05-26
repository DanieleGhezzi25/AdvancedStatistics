library(cmdstanr) # we use STAN
library(coda)
library(posterior)

set.seed(42)
data <- read.csv("ex3/logistic_data.csv")
summary(data)
data <- list(
  N = nrow(data),
  x = data$x,
  y = data$y
)
sprintf("Number of observations: %d", data$N)


# 1. Posterior formulation and implementation
mod <- cmdstan_model("ex3/logistic_regression.stan")


# 2. MCMC fitting and diagnostics
fit <- mod$sample(data=data,
    chains=4,
    parallel_chains=4,
    iter_warmup=5000,
    iter_sampling=50000,
    refresh = 0, # Silences the iteration progress updates
)
fit$summary(c("beta0", "beta1"))
fit$cmdstan_diagnose()

mcmc_full <- cmdstanr::as_mcmc.list(fit)
mcmc_chain <- mcmc_full[, c("beta0", "beta1")]
summary(mcmc_chain)

# autocorrelation
mcmc_mat <- as.matrix(mcmc_chain)
par(mfrow=c(1, 2)) 
for (j in 1:ncol(mcmc_mat)) {
    acf(mcmc_mat[,j], 
        main=paste("ACF for", colnames(mcmc_mat)[j]))
}

# ESS
ess <- effectiveSize(mcmc_mat)
print(ess)

# trace plots
par(mfrow=c(1,2))
for (j in 1:ncol(mcmc_mat)) {
    plot(mcmc_mat[,j], type="l",
    main=paste("Trace theta", j),
    ylab="", xlab="")
}
