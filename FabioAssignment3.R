library(cmdstanr)
library(posterior)
library(ggplot2)
library(bayesplot)

data <- read.csv('logistic_data.csv')
x <- data$x
y <- data$y
N <- length(x)

x_min <- min(x)
x_max <- max(x)
x_mean <- mean(x)
x_sd <- sd(x)
n_success <- sum(y)
n_failure <- N - n_success

cat("Number of observations:", N, "\n")
cat("x: min =", x_min, ", max =", x_max, ", mean =", x_mean, ", sd =", x_sd, "\n")
cat("y: successes =", n_success, ", failures =", n_failure, "\n")

pdf("data_scatter.pdf", width = 6, height = 4)
ggplot(data, aes(x = x, y = y)) +
  geom_point(alpha = 0.6) +
  labs(title = "Scatter plot of y against x",
       x = "x",
       y = "y") +
  theme_minimal()
dev.off()

# Histogram of x
pdf("x_histogram.pdf", width = 6, height = 4)
ggplot(data, aes(x = x)) +
  geom_histogram(bins = 20) +
  labs(title = "Histogram of x", x = "x", y = "Count") +
  theme_minimal()
dev.off()



stan_data <- list(N=N, x=x, y=y)
model1 <- cmdstan_model("model1.stan")
fit1 <- model1$sample(data=stan_data, chains=4, parallel_chains=4, iter_warmup=5000, iter_sampling=2000, seed=42)

fit1$summary(c("beta0", "beta1"))
fit1$cmdstan_diagnose()

# Draw posterior samples
draws <- fit1$draws(c("beta0", "beta1"))
draws_df <- as_draws_df(draws)
beta0_draws = draws_df$beta0
beta1_draws = draws_df$beta1

draws_summary <- posterior::summarise_draws(draws, mean, median, sd, ~posterior::quantile2(.x, probs = c(0.025, 0.975)))
print(draws_summary)
b0_draws_summary <- subset(draws_summary, variable == "beta0")
b1_draws_summary <- subset(draws_summary, variable == "beta1")

# Trace plots
pdf("trace_plot1.pdf", width = 8, height = 3)
mcmc_trace(draws, pars = c("beta0"))
dev.off()
pdf("trace_plot2.pdf", width = 8, height = 3)
mcmc_trace(draws, pars = c("beta1"))
dev.off()

# ACF plots
pdf("acf_plots.pdf", width = 8, height = 6)
mcmc_acf(draws, pars = c("beta0", "beta1"))
dev.off()


# Marginal posterior distributions
pdf("marginal_posterior1.pdf", width = 4, height = 4)
ggplot(draws_df, aes(x=beta0)) +
  geom_histogram(aes(y=after_stat(density)), bins=50) +
  geom_density() +
  geom_vline(xintercept = b0_draws_summary$mean, linetype = "solid") +
  geom_vline(xintercept = b0_draws_summary$q2.5, linetype = "dashed") +
  geom_vline(xintercept = b0_draws_summary$q97.5, linetype = "dashed") +
  labs(title = expression("Posterior distribution of" ~ beta[0]),
       x = expression(beta[0]),
       y = "Density") +
  theme_minimal()
dev.off()

pdf("marginal_posterior2.pdf", width = 4, height = 4)
ggplot(draws_df, aes(x=beta1)) +
  geom_histogram(aes(y=after_stat(density)), bins=50) +
  geom_density() +
  geom_vline(xintercept = b1_draws_summary$mean, linetype = "solid") +
  geom_vline(xintercept = b1_draws_summary$q2.5, linetype = "dashed") +
  geom_vline(xintercept = b1_draws_summary$q97.5, linetype = "dashed") +
  labs(title = expression("Posterior distribution of" ~ beta[1]),
       x = expression(beta[1]),
       y = "Density") +
  theme_minimal()
dev.off()


# Joint posteriors
vars <- c("beta0", "beta1")
cor(as.matrix(draws_df[, vars]))

pdf("joint_density.pdf", width = 5, height = 4)
ggplot(draws_df, aes(x = beta0, y = beta1)) + geom_density_2d_filled() + theme_minimal() + labs(
  x = expression(beta[0]), y = expression(beta[1]))
dev.off()

pdf("joint_scatter.pdf", width = 4, height = 4)
ggplot(draws_df, aes(x=beta0, y=beta1)) + geom_point(alpha=0.3) + geom_density_2d() + labs(
  x = expression(beta[0]), y = expression(beta[1]))
dev.off()




# New covariate value
x_new <- 2.0

pi_new_draws <- plogis(beta0_draws + beta1_draws * x_new)
pi_new_mean <- mean(pi_new_draws)
pi_new_median <- median(pi_new_draws)
pi_new_sd <- sd(pi_new_draws)
pi_new_ci <- quantile(pi_new_draws, probs = c(0.025, 0.975))

# Posterior predictive probability
pred_prob <- mean(pi_new_draws)

cat("Posterior mean:", pi_new_mean, "\n")
cat("Posterior median:", pi_new_median, "\n")
cat("Posterior SD:", pi_new_sd, "\n")
cat("95% credible interval: (", pi_new_ci[1], ", ", pi_new_ci[2], ")\n")
cat("Posterior predictive probability P(y_new=1 | x_new,D):", pred_prob, "\n")

pdf("pi_new_histogram.pdf", width = 7, height = 5)
pi_new_df <- data.frame(pi_new = pi_new_draws)
ggplot(pi_new_df, aes(x = pi_new)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30) +
  geom_density() +
  geom_vline(xintercept = pi_new_mean, linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = pi_new_ci[1], linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = pi_new_ci[2], linetype = "dashed", linewidth = 1) +
  labs(
    title = expression("Posterior distribution of " ~ pi[new]),
    x = expression(pi[new]),
    y = "Density"
  ) +
  theme_minimal()
dev.off()




#https://mc-stan.org/loo/articles/online-only/faq.html

library(loo)

model2 <- cmdstan_model("model2.stan")

fit2 <- model2$sample(data=stan_data, chains=4, parallel_chains=4, iter_warmup=1000, iter_sampling=2000, seed=42)

fit2$summary(c("beta0"))
fit2$cmdstan_diagnose()

draws2 <- fit2$draws(c("beta0"))
mcmc_trace(draws2, pars = c("beta0"))
mcmc_acf(draws2, pars = c("beta0"))


log_lik1 <- fit1$draws("log_lik")
log_lik2 <- fit2$draws("log_lik")

loo1 <- loo(log_lik1)
loo2 <- loo(log_lik2)

loo_compare(loo1, loo2)
