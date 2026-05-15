# load data
filepath = '/home/dghezzi/AdvancedStatistics/ex2/assignment2_changepoint_counts.csv'
data = read.csv(filepath)
set.seed(4) # for reproducibility

# summary statistics
summary(data$x)
sprintf("Variance: %.2f", var(data$x))

alpha_prior <- 1
beta_prior <- 0.1

# 2
# 2.a write the three update steps of one Gibbs iteration
gibbs_sampling_posterior <- function(n_steps, data, alpha, beta, start = c(1, 1, 10)) {
    
    # initialize empty vectors
    lambda1 <- numeric(n_steps)
    lambda2 <- numeric(n_steps)
    m <- numeric(n_steps)
    
    # initialize first values
    lambda1[1] <- start[1]
    lambda2[1] <- start[2]
    m[1] <- start[3]
    
    sumN <- sum(data$x)
    N <- length(data$x)
    Mset <- seq(5, N-5, 1)
    
    for (t in 2:n_steps) {
        current_m <- m[t-1]
        sumM_current <- sum(data$x[1:current_m])
        
        # lambda updates
        lambda1[t] <- rgamma(1, shape = alpha + sumM_current, rate = beta + current_m)
        lambda2[t] <- rgamma(1, shape = alpha + sumN - sumM_current, rate = beta + N - current_m)
        
        # m update
        logm <- numeric(length(Mset))
        
        for (i in 1:length(Mset)) {
        m_cand <- Mset[i] 
        sumM_cand <- sum(data$x[1:m_cand]) 
        
        logm[i] <- sumM_cand * log(lambda1[t]) - m_cand * lambda1[t] + (sumN - sumM_cand) * log(lambda2[t]) - (N - m_cand) * lambda2[t]
        }
        
        weights <- exp(logm - max(logm))
        probs_m <- weights / sum(weights)
        
        m[t] <- sample(x = Mset, size = 1, prob = probs_m)
  }
  
  print("Gibbs sampling completed.")

  return(cbind(lambda1 = lambda1, lambda2 = lambda2, m = m))
}

# 2.b Run three chains starting from different values and show the trace plots 

starts <- list(
  c(4, 5, 10), 
  c(6, 4, 20), 
  c(2, 6, 25)
)

traces <- vector("list", length(starts))
n_steps <- 50000

for (i in 1:length(starts)) {
  traces[[i]] <- gibbs_sampling_posterior(n_steps = n_steps, 
                                 data = data, 
                                 alpha = alpha_prior,
                                 beta = beta_prior,
                                 start = starts[[i]])
}

par(mfrow = c(3, 3), mar = c(2.5, 2.5, 2, 1), mgp = c(1.5, 0.5, 0))
chain_cols <- c("blue", "red", "forestgreen")
param_names <- c("lambda1", "lambda2", "m")
param_exprs <- c(expression(lambda[1]), expression(lambda[2]), expression(m))

for (p in 1:3) {
  
  param <- param_names[p]
  y_label <- param_exprs[[p]]
  
  global_ylim <- range(c(traces[[1]][, param], traces[[2]][, param], traces[[3]][, param]))
  
  for (i in 1:3) { 
    
    plot_title <- paste(param_names[p], "- Chain", i)
    
    plot(traces[[i]][, param], type = "l", col = chain_cols[i], 
         ylim = global_ylim, main = plot_title, 
         xlab = "Iteration", ylab = y_label,
         cex.main = 1.2, cex.lab = 1) 
  }
}
par(mfrow = c(1, 1))

# 2.c burn-in period?
burn_in <- 5000
sample_gibbs <- vector("list", length(traces))
for (i in 1:length(traces)) {
  sample_gibbs[[i]] <- traces[[i]][(burn_in + 1):n_steps, ]
}

# plot with and without burn-in ? i think it is not possible since we are in three dimensions? maybe scatterplot in 3d to see the burn-in

# plot histogram of the samples
lambda1_samples <- sample_gibbs[[1]][, "lambda1"]
lambda2_samples <- sample_gibbs[[1]][, "lambda2"]
m_samples <- sample_gibbs[[1]][, "m"]
par(mfrow = c(1, 3), mar = c(2.5, 2.5, 3.5, 1))
hist(lambda1_samples, breaks = 30, main = expression(lambda[1]), xlab = expression(lambda[1]), col = "lightblue", border = "white")
hist(lambda2_samples, breaks = 30, main = expression(lambda[2]), xlab = expression(lambda[2]), col = "lightblue", border = "white")
hist(m_samples, breaks = 30, main = "m", xlab = "m", col = "lightblue", border = "white")
par(mfrow = c(1, 1))


# 2.d Autocorrelation? Effective sample size?
lagMAX <- 100
par(mfrow = c(1, 3), mar = c(2.5, 2.5, 3.5, 1))
acf_lambda1 <-acf(lambda1_samples, lag.max = lagMAX, main = "ACF for lambda1")
acf_lambda2 <- acf(lambda2_samples, lag.max = lagMAX, main = "ACF for lambda2")
acf_m <- acf(m_samples, lag.max = lagMAX, main = "ACF for m")
par(mfrow = c(1, 1))

calculate_tau_int <- function(acf) {

  lags <- as.numeric(acf$lag)
  rho <- as.numeric(acf$acf)
  positive_lags <- which(rho[-1] > 0)
  
  if (length(positive_lags) == 0) {
    return(1)
  } else {
    last_positive <- positive_lags[1]
    for (i in seq_along(positive_lags)) {
      if (positive_lags[i] == i) {
        last_positive <- positive_lags[i]
      } else {
        break
      }
    }
    return(1 + 2 * sum(rho[2:(last_positive + 1)]))
  }
}

tau_int_lambda1 <- calculate_tau_int(acf_lambda1)
tau_int_lambda2 <- calculate_tau_int(acf_lambda2)
tau_int_m <- calculate_tau_int(acf_m)

N_samples <- length(lambda1_samples)
ess_lambda1 <- N_samples / (2*tau_int_lambda1)
ess_lambda2 <- N_samples / (2*tau_int_lambda2)
ess_m <- N_samples / (2*tau_int_m)
print(paste("ESS for lambda1:", round(ess_lambda1, 0)))
print(paste("ESS for lambda2:", round(ess_lambda2, 0)))
print(paste("ESS for m:", round(ess_m, 0)))


# 2.e For each par: posterior mean, posterior median, posterior standard deviation, MAP estimate, and a 95% credible interval?
get_hpd_interval <- function(samples, prob = 0.95) {

  # the HPD interval is the narrowest interval that contains 95% of the posterior mass
  
  sorted_samples <- sort(samples)
  n <- length(sorted_samples)
  
  # number of samples needed to cover the target probability
  window_size <- round(n * prob)
  
  # define the start and end indices for all possible sliding windows across the sorted samples
  start_idx <- 1:(n - window_size)
  end_idx <- start_idx + window_size
  
  # calculate the physical width of every single window
  window_widths <- sorted_samples[end_idx] - sorted_samples[start_idx]
  
  # find the index of the narrowest window
  best_idx <- which.min(window_widths)
  
  # Return the boundaries of that narrowest window
  return(c(lower = sorted_samples[best_idx], upper = sorted_samples[end_idx[best_idx]]))
}

posterior_summary <- function(samples) {
  mean_val <- mean(samples)
  median_val <- median(samples)
  sd_val <- sd(samples)
  ci <- get_hpd_interval(samples, prob = 0.95)
  return(list(mean = mean_val, median = median_val, sd = sd_val, ci = ci))
}

lambda1_summary <- posterior_summary(lambda1_samples)
lambda2_summary <- posterior_summary(lambda2_samples)
m_summary <- posterior_summary(m_samples)

print("Lambda1 Summary:")
print(lambda1_summary)
print("Lambda2 Summary:")
print(lambda2_summary)
print("m Summary:")
print(m_summary)

# we are calculating central credible intervals, not HPD intervals
# to get the MAP estimate, there are two options:
# 1) we calculate the log-posterior for each sample and we take the one with the highest value;
# 2) we use KDE to estimate the density and we find its maximum. 
# method 1
calculate_log_posterior <- function(lambda1, lambda2, m, data, alpha, beta) {
  N <- length(data$x)
  
  sumN <- sum(data$x)
  sumM <- sum(data$x[1:m])
  
  log_prior_lambda1 <- dgamma(lambda1, shape = alpha, rate = beta, log = TRUE)
  log_prior_lambda2 <- dgamma(lambda2, shape = alpha, rate = beta, log = TRUE)
  
  log_likelihood <- sumM * log(lambda1) - m * lambda1 + (sumN - sumM) * log(lambda2) - (N - m) * lambda2
  
  return(log_prior_lambda1 + log_prior_lambda2 + log_likelihood)
}

log_posteriors <- mapply(calculate_log_posterior, 
                         lambda1 = lambda1_samples, 
                         lambda2 = lambda2_samples, 
                         m = m_samples, 
                         MoreArgs = list(data = data, alpha = alpha_prior, beta = beta_prior))
max_index <- which.max(log_posteriors)
map_lambda1 <- lambda1_samples[max_index]
map_lambda2 <- lambda2_samples[max_index]
map_m <- m_samples[max_index]
sprintf("MAP estimate for lambda1: %.4f", map_lambda1)
sprintf("MAP estimate for lambda2: %.4f", map_lambda2)
sprintf("MAP estimate for m: %.0f", map_m)

# method 2 (no KDE on m since it is discrete)
kde_lambda1 <- density(lambda1_samples)
kde_lambda2 <- density(lambda2_samples)
map_lambda1_kde <- kde_lambda1$x[which.max(kde_lambda1$y)]
map_lambda2_kde <- kde_lambda2$x[which.max(kde_lambda2$y)]
sprintf("MAP estimate for lambda1 (KDE): %.4f", map_lambda1_kde)
sprintf("MAP estimate for lambda2 (KDE): %.4f", map_lambda2_kde)

# plot kde and map estimates with histogram of the samples
par(mfrow = c(1, 3), mar = c(2.5, 2.5, 3.5, 1))

hist(lambda1_samples, breaks = 30, main = expression(lambda[1]), xlab = expression(lambda[1]), col = "lightblue", border = "white", freq = FALSE)
lines(kde_lambda1, col = "red", lwd = 2)
abline(v = map_lambda1, col = "blue", lwd = 2, lty = 2)
abline(v = map_lambda1_kde, col = "#1e9600", lwd = 2, lty = 3)
abline(v = lambda1_summary$ci[1], col = "purple", lwd = 2, lty = 4)
abline(v = lambda1_summary$ci[2], col = "purple", lwd = 2, lty = 4)
legend("topright", legend = c("MAP KDE", "MAP samples", "95% CI (HPD)"), col = c("#1e9600", "blue", "purple"), lwd = 2, lty = c(3, 2, 4), cex = 0.8)

hist(lambda2_samples, breaks = 30, main = expression(lambda[2]), xlab = expression(lambda[2]), col = "lightblue", border = "white", freq = FALSE)
lines(kde_lambda2, col = "red", lwd = 2)
abline(v = map_lambda2, col = "blue", lwd = 2, lty = 2)
abline(v = map_lambda2_kde, col = "#1e9600", lwd = 2, lty = 3)
abline(v = lambda2_summary$ci[1], col = "purple", lwd = 2, lty = 4)
abline(v = lambda2_summary$ci[2], col = "purple", lwd = 2, lty = 4)
legend("topright", legend = c("MAP KDE", "MAP samples", "95% CI (HPD)"), col = c("#1e9600", "blue", "purple"), lwd = 2, lty = c(3, 2, 4), cex = 0.8)

hist(m_samples, breaks = 30, main = "m", xlab = "m", col = "lightblue", border = "white", freq = FALSE)
abline(v = map_m, col = "blue", lwd = 2, lty = 2)
abline(v = m_summary$ci[1], col = "purple", lwd = 2, lty = 4)
abline(v = m_summary$ci[2], col = "purple", lwd = 2, lty = 4)
legend("topright", legend = c("MAP samples", "95% CI (HPD)"), col = c("blue", "purple"), lwd = 2, lty = c(2, 4), cex = 0.8)

par(mfrow = c(1, 1))

# 2.f Estimate the posterior probability P(lambda2 > lambda1 | D)
# approx of double integral over the joint posterior of lambda1 and lambda2 where lambda2 > lambda1 as integration domain
prob_lambda2_greater <- mean(lambda2_samples > lambda1_samples)
print(paste("P(lambda2 > lambda1 | D) =", round(prob_lambda2_greater, 4)))


# 3.b Estimate the posterior predictive mean
n_predictive_samples <- 1000 * length(lambda1_samples)
sample_predictive_distribution <- rpois(n = n_predictive_samples, lambda = lambda1_samples)
print(paste("Posterior predictive mean:", round(mean(sample_predictive_distribution), 4)))


# 3.c Estimate the probability that the next count is greater than twice the mean of the observed counts
sample_mean <- mean(data$x)
prob_twice_mean <- mean(sample_predictive_distribution > 2*sample_mean)
print(paste("Probability that the next count is greater than twice the mean of the observed counts:", round(prob_twice_mean, 4)))

# 3.d Histogram of posterior predictive distribution, with sample_mean and 2*sample_mean
par(mfrow = c(1, 1))
breaks = seq(min(sample_predictive_distribution) - 0.5, 
             max(sample_predictive_distribution) + 0.5, 
             by = 1) # ensure that each integer count gets its own bin
hist(sample_predictive_distribution, breaks = breaks, main = "Posterior Predictive Distribution", xlab = "Predicted Count", col = "lightblue", border = "white", freq = FALSE)
abline(v = sample_mean, col = "red", lwd = 2, lty = 2)
abline(v = 2*sample_mean, col = "green", lwd = 2, lty = 2)
legend("topright", legend = c("Mean of observed counts", "Twice the mean of observed counts"), col = c("red", "green"), lwd = 2, lty = 2, cex = 0.8)
par(mfrow = c(1, 1))


# 4. Model comparison
N <- length(data$x)
sumN <- sum(data$x)
Mset <- seq(5, N-5, 1)
len_M <- length(Mset)

# 4.a Marginal likelihood for M0 (on log scale, ignoring the 1/prod(x_i!) term)
log_marg_lik_M0 <- alpha_prior * log(beta_prior) - 
                   lgamma(alpha_prior) + 
                   lgamma(alpha_prior + sumN) - 
                   (alpha_prior + sumN) * log(beta_prior + N)
print(log_marg_lik_M0)

# 4.b Marginal likelihood for M1 (on log scale, ignoring the 1/prod(x_i!) term)
log_terms_M1 <- numeric(len_M)

for (i in 1:len_M) {
  m <- Mset[i]
  sumM <- sum(data$x[1:m])
  log_terms_M1[i] <- 2 * alpha_prior * log(beta_prior) - 
                     log(len_M) - 
                     2 * lgamma(alpha_prior) + 
                     lgamma(alpha_prior + sumM) + 
                     lgamma(alpha_prior + sumN - sumM) - 
                     (alpha_prior + sumM) * log(beta_prior + m) - 
                     (alpha_prior + sumN - sumM) * log(beta_prior + N - m)
}

log_marg_lik_M1 <- log(sum(exp(log_terms_M1)))
print(log_marg_lik_M1)

# 4.c Compute the log Bayes factor
log_BF10 <- log_marg_lik_M1 - log_marg_lik_M0
sprintf("Log Bayes Factor (log BF10): %.4f", log_BF10)
