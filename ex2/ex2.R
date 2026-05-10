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
  c(1, 1, 10), 
  c(5, 5, 20), 
  c(0.1, 0.1, 30)
)

traces <- vector("list", length(starts))
n_steps <- 50000

for (i in 1:length(starts)) {
  traces[[i]] <- gibbs_posterior(n_steps = n_steps, 
                                 data = data, 
                                 alpha = alpha_prior,
                                 beta = beta_prior,
                                 start = starts[[i]])
}