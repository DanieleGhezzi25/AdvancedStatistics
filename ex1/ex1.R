library(ggplot2)

# Load the data into data frame
file <- "/home/dghezzi/AdvancedStatistics/ex1/dataset.txt"
data <- read.table(file, header = FALSE)
colnames(data) <- c("data")

# 1.a Derivation of posterior distribution
# Since as prior we have a Gamma distribution, the posterior distribution will also be a Gamma distribution.
# The parameters of the posterior distribution can be calculated as follows:

alpha_prior <- 0.001
beta_prior <- 0.001

alpha_post <- alpha_prior + sum(data$data)
beta_post <- beta_prior + nrow(data)
post_mean <- alpha_post / beta_post # mean of the posterior distribution
post_var <- alpha_post / (beta_post^2) # variance of the posterior distribution


cat("The posterior distribution is Gamma with alpha =", alpha_post, "and beta =", beta_post, "\n")
cat("The mean of the posterior distribution is:", post_mean, "\n")
cat("The variance of the posterior distribution is:", post_var, "\n")


# 1.b Maximum a posteriori (MAP) estimate of lambda
# From analytical calculations, the MAP estimate of lambda can be calculated as follows:

map_estimate <- (alpha_post - 1) / beta_post
cat("The MAP estimate of lambda is:", map_estimate, "\n")

# Plot of the posterior distribution and the MAP estimate
x <- seq(0, 10, length.out = 1000)
posterior <- dgamma(x, shape = alpha_post, rate = beta_post)
plot_data <- data.frame(lambda = x, density = posterior)
ggplot(plot_data, aes(x = lambda, y = density)) +
  geom_area(fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_vline(xintercept = map_estimate, color = "red", linetype = "dashed", linewidth = 0.5) +
  labs(title = "Posterior Distribution",
       subtitle = paste("MAP Estimate of lambda =", round(map_estimate, 2)),
       x = expression(lambda), 
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
ggsave("ex1/posterior_MAPestimate_Poisson.png", width = 8, height = 6, dpi = 300)

# 1.c Provide an adequate estimate of the uncertainty on lambda
# We check skewness and kurtosis of the posterior distribution to see if a gaussian approximation is reasonable
skewness <- 2 / sqrt(alpha_post)
kurtosis <- 6 / alpha_post
cat("The skewness of the posterior distribution is:", skewness, "\n")
cat("The kurtosis of the posterior distribution is:", kurtosis, "\n")

# As skewness and kurtosis are close to 0, we can use a Gaussian approximation to estimate the uncertainty on lambda
var_lambda <- (alpha_post - 1) / (beta_post^2)
cat("The variance of the MAP estimate of lambda is:", var_lambda, "\n")
cat("The MAP estimate of lambda with gaussian approximation is:", map_estimate, "±", sqrt(var_lambda), "\n")


# 2.a Posterior predictive distribution for a new observation
# The posterior predictive distribution is a Negative Binomial pmf with parameters r = alpha_post and p = beta_post / (beta_post + 1)
# where the k successes correspond to x_new

prob_param <- beta_post / (beta_post + 1)
x_max <- qnbinom(0.999, size = alpha_post, prob = prob_param) # sensible maximum
x_vals <- 0:x_max # sequence of integers as we are dealing with a discrete distribution
pmf <- dnbinom(x_vals, size = alpha_post, prob = prob_param)

plot_data <- data.frame(
  x_new = x_vals, 
  probability = pmf
)

# plot
ggplot(plot_data, aes(x = x_new, y = probability)) +
  geom_segment(aes(x = x_new, xend = x_new, y = 0, yend = probability), color = "orange", linewidth = 3) +
  labs(title = "Posterior Predictive Distribution",
       x = expression(x[new] ~ "(Number of Events)"),
       y = "Probability Mass") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
ggsave("ex1/posteriorPredictive_Poisson.png", width = 8, height = 6, dpi = 300)

# 2.b Compute the probability that the next observation exceeds twice the empirical average
threshold <- 2 * mean(data$data)
exceeds_prob <- 1 - pnbinom(q = threshold, size = alpha_post, prob = beta_post / (beta_post + 1))
cat("The probability that the next observation exceeds twice the empirical average is:", exceeds_prob, "\n")


# 3.a Consider an alternative model with a Negative Binomial distribution for the data, with parameters size r and prob pi
# Give the posterior distribution for the parameters by choosing appropriate priors

# We assume a Gamma prior for the size parameter and a Beta prior for the prob parameter. The likelihood is given by the Negative Binomial distribution.
# The posterior distribution can be calculated using Bayes' theorem, but it does not have a closed-form solution due to the complexity of the evidence integral.
# However, we can still plot the unnormalized posterior to visualize the distribution of the parameters.
# As probabilities could be infinitesimal, we will plot the log of the unnormalized posterior distribution.

# we define negative log-posterior function as later in point 3.b we will use it for optimization
neg_log_post <- function(params, data) {
  r <- params[1]
  pi <- params[2]
  
  # Log-Likelihood of the Negative Binomial distribution
  log_lik <- sum(dnbinom(data, size = r, prob = pi, log = TRUE))
  
  # Weakly informative Gamma for r (shape=0.01, rate=0.01)
  log_prior_r <- dgamma(r, shape = 0.01, rate = 0.01, log = TRUE)
  
  # Jeffreys prior log-Beta for pi (shape1=1/2, shape2=1/2)
  log_prior_pi <- dbeta(pi, shape1 = 1/2, shape2 = 1/2, log = TRUE)

  log_post <- log_lik + log_prior_r + log_prior_pi
  return(-log_post)
}

# We set up the grid for plotting
sample_mean <- mean(data$data)
sample_var <- var(data$data)

# This estimated values help us to set a grid that is focused where the posterior is signifcant
pi_mom <- sample_mean / sample_var
r_mom <- (sample_mean^2) / (sample_var - sample_mean)

r_vals <- seq(max(0.01, r_mom * 0.1), r_mom * 2, length.out = 300) # we use the max as r cannot be zero
pi_vals <- seq(max(0.01, pi_mom - 0.2), min(0.99, pi_mom + 0.2), length.out = 300)
grid <- expand.grid(r = r_vals, pi = pi_vals)

grid$log_post <- apply(grid, 1, function(row) {
  -neg_log_post(c(row["r"], row["pi"]), data = data$data)
})

# Shift maximum to 0 to prevent underflow, then exponentiate to recover the posterior values
grid$posterior <- exp(grid$log_post - max(grid$log_post))

ggplot(grid, aes(x = r, y = pi, z = posterior)) +
  geom_contour_filled(bins = 20) +
  labs(title = "2D Unnormalized Posterior Distribution",
       x = expression("Size Parameter " ~ r),
       y = expression("Probability Parameter " ~ pi),
       fill = "Relative\nDensity") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")
ggsave("ex1/posterior_MAPestimate_NegativeBinomial.png", width = 8, height = 6, dpi = 300)

# 3.2 Assuming gaussian shape, compute MAP estimates and uncertainties for the parameters
# We use numerical optimization to find the MAP estimates and the Hessian to estimate the covariance matrix

# We use "L-BFGS-B" method as we want boundaries: r > 0 and 0 < pi < 1
initial_params <- c(r_mom, pi_mom)
fit <- optim(par = initial_params, 
             fn = neg_log_post, 
             data = data$data, 
             method = "L-BFGS-B", 
             lower = c(1e-5, 1e-5), # we set lower bounds to (r, pi)
             upper = c(Inf, 1 - 1e-5), # we set upper bounds to (r, pi)
             hessian = TRUE)

map_r <- fit$par[1]
map_pi <- fit$par[2]

# find the covariance matrix by inverting the Hessian matrix
cov_matrix <- solve(fit$hessian) 

std_r <- sqrt(cov_matrix[1, 1])
std_pi <- sqrt(cov_matrix[2, 2])

cat("MAP Estimate for r:", map_r, "±", std_r, "\n")
cat("MAP Estimate for pi:", map_pi, "±", std_pi, "\n")
cat("\nCovariance Matrix:\n")
print(cov_matrix)


# 3.3 Evidence in Gaussian approximation

log_evidence_nb <- -fit$value + log(2 * pi) + 0.5 * log(det(cov_matrix))
cat("The log-evidence in Gaussian approximation is:", log_evidence_nb, "\n")
# cat("The evidence in Gaussian approximation is:", exp(log_evidence), "\n") => underflow!


# 3.4 Calculate the Bayes factor between the two models
log_evidence_Poisson <- alpha_prior * log(beta_prior) + 
                        lgamma(alpha_post) - 
                        lgamma(alpha_prior) - 
                        alpha_post * log(beta_post) - 
                        sum(lfactorial(data$data))

BF <- exp(log_evidence_nb - log_evidence_Poisson)
cat("The Bayes factor in favor of the Negative Binomial model over the Poisson model is:", BF, "\n")
