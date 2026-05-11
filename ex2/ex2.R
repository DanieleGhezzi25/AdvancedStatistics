# load data
filepath = '/home/dghezzi/AdvancedStatistics/ex2/assignment2_changepoint_counts.csv'
data = read.csv(filepath)

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

# 2.d Autocorrelation? Effective sample size?
par(mfrow = c(1, 3), mar = c(2.5, 2.5, 3.5, 1))
acf(traces[[1]][, "lambda1"], lag.max = 100, main = "ACF for lambda1")
acf(traces[[1]][, "lambda2"], lag.max = 100, main = "ACF for lambda2")
acf(traces[[1]][, "m"], lag.max = 100, main = "ACF for m")
par(mfrow = c(1, 1))

lagMAX <- 15 # from the plot we see that the chain is not autocorrelated starting form lag=15. We remove the subsequent lag since it is noise (it should be zero)
acf_values_lambda1 <- acf(sample_gibbs[[1]][, "lambda1"], plot = FALSE, lag.max = lagMAX)$acf
acf_values_lambda2 <- acf(sample_gibbs[[1]][, "lambda2"], plot = FALSE, lag.max = lagMAX)$acf
acf_values_m <- acf(sample_gibbs[[1]][, "m"], plot = FALSE, lag.max = lagMAX)$acf

tau_int_lambda1 <- 1 + 2 * sum(acf_values_lambda1[-1]) # we remove the first autocorrelation as it is 1
tau_int_lambda2 <- 1 + 2 * sum(acf_values_lambda2[-1])
tau_int_m <- 1 + 2 * sum(acf_values_m[-1])

N_samples <- nrow(sample_gibbs[[1]])
ess_lambda1 <- N_samples / (2*tau_int_lambda1)
ess_lambda2 <- N_samples / (2*tau_int_lambda2)
ess_m <- N_samples / (2*tau_int_m)
print(paste("ESS for lambda1:", round(ess_lambda1, 0)))
print(paste("ESS for lambda2:", round(ess_lambda2, 0)))
print(paste("ESS for m:", round(ess_m, 0)))


# 2.e For each par: posterior mean, posterior median, posterior standard deviation, MAP estimate, and a 95% credible interval?
lambda1_samples <- sample_gibbs[[1]][, "lambda1"]
lambda2_samples <- sample_gibbs[[1]][, "lambda2"]
m_samples <- sample_gibbs[[1]][, "m"]
posterior_summary <- function(samples) {
  mean_val <- mean(samples)
  median_val <- median(samples)
  sd_val <- sd(samples)
  ci_lower <- quantile(samples, 0.025) 
  ci_upper <- quantile(samples, 0.975)
  
  return(list(mean = mean_val, median = median_val, sd = sd_val, ci = c(ci_lower, ci_upper)))
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
# 2) we use KDE to estimate the density and then we run an optimization to find its maximum. 