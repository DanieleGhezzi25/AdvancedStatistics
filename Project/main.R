# ================================================================ #
#                        PREPROCESSING                             #
# ================================================================ #
library(dplyr)
library(ggplot2)
library(cmdstanr)
library(posterior)
library(bayesplot)

path <- "co2_mm_mlo.txt"
df <- read.table(path, comment.char = "#", col.names = c(
  "year",
  "month",
  "decimal_date",
  "monthly_average",
  "de_seasonalized",
  "ndays",
  "stdev",
  "unc_monthly_mean"
))

# Monthly CO2 plot
monthly_plot <- ggplot(df, aes(decimal_date, monthly_average)) +
  geom_line(linewidth = 0.3) +
  geom_point(size = 0.5) +
  labs(
    x = "Year",
    y = expression(CO[2] ~ "(ppm)"),
    title = "Monthly CO2 measurements"
  ) +
  theme_minimal()

ggsave("monthly_co2.pdf", monthly_plot,width = 8,height = 5)

# filter problematic months
df <- df %>%
  filter(
    !(year == 1984 & month == 4),
    !(year == 1975 & month == 12)
  )

# Annual average
annual <- df %>%
  group_by(year) %>%
  summarise(
    annual_mean =
      if (first(year) <= 1974) {
        mean(monthly_average)
      } else {
        weighted.mean(monthly_average, ndays)
      },
    
    annual_unc =
      if (first(year) <= 1974) {
        NA
      } else {
        sqrt(sum(ndays^2 * unc_monthly_mean^2)) / sum(ndays)
      },
    
    n_months = sum(!is.na(monthly_average)),
    
    .groups = "drop"
  )


annual_plot <- ggplot(annual, aes(year, annual_mean)) +
  geom_line() +
  geom_point() +
  geom_errorbar(
    data = subset(annual, !is.na(annual_unc)),
    aes(
      ymin = annual_mean - annual_unc,
      ymax = annual_mean + annual_unc
    ),
    width = 0.3
  ) +
  labs(
    x = "Year",
    y = expression(CO[2] ~ "(ppm)"),
    title = "Annual average CO2"
  ) +
  theme_minimal()

ggsave("annual_co2.pdf", annual_plot, width = 8, height = 5)

print(as.data.frame(annual))
print(min(annual$annual_mean))
print(max(annual$annual_mean))
print(min(annual$year))
print(max(annual$year)-min(annual$year))



# ================================================================ #
#                      LINEAR MODEL FITTING                        #
# ================================================================ #

library(cmdstanr)
library(bayesplot)
library(GGally)

analyze_MCMC <- function(fit, model, stan_data, param_names, fit_name) {
  # NOTE: param_names must leave sigma as the last element
  # Compute statistics and plots for the specified MCMC fit and parameters
  
  dir.create(fit_name, showWarnings = FALSE)
  
  draws <- fit$draws(param_names)
  draws_df <- as_draws_df(draws)
  
  # 1. Compute summary statistics
  summary_stats <- posterior::summarise_draws(
    draws,
    mean,
    median,
    sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975)),
    rhat,
    ess_bulk,
    ess_tail
  )
  print(summary_stats)
  
  # 2. MAP estimate
  init_fun <- function() {
    inits <- list(
      b0 = as.numeric(fit$summary("b0")$mean),
      alpha = as.numeric(fit$summary("alpha")$mean),
      logsigma = as.numeric(fit$summary("logsigma")$mean) 
    )
    # if the model is quadratic, then:
    if ("b2" %in% fit$summary()$variable) {
      inits$b2 <- as.numeric(fit$summary("b2")$mean)
    }
    return(inits)
  }
  
  map_fit <- model$optimize(
    data = stan_data, 
    seed = 42, 
    init = init_fun
  )
  
  cat("\n=== MAP estimates & 95% CI ===\n")
  for (param in param_names) {
    map_val <- map_fit$summary(param)$estimate
    q2.5 <- summary_stats$q2.5[summary_stats$variable == param]
    q97.5 <- summary_stats$q97.5[summary_stats$variable == param]
    
    cat(sprintf("%s: %.3f  (95%% CI: [%.3f, %.3f])\n", 
                param, map_val, q2.5, q97.5))
  }
  
  # 3. Trace plots
  for (i in seq_along(param_names)) {
    param <- param_names[i]
    pdf(file.path(fit_name, paste0("trace_plot_", fit_name, "_", param, ".pdf")), width = 8, height = 3)
    print(mcmc_trace(draws, pars = param))
    dev.off()
  }
  
  # 4. ACF plots
  pdf(file.path(fit_name, paste0("acf_plots_", fit_name, ".pdf")), width = 8, height = 6)
  print(mcmc_acf(draws, pars = param_names))
  dev.off()
  
  
  # 6. Joint posteriors depending on number of parameters
  corr_mat <- cor(as.matrix(draws_df[, param_names]))
  print(corr_mat)
  #pairs <- combn(param_names, 2, simplify = FALSE) # combination of all variables taken 2 at a time
  #for (pair in pairs) {
  #  var1 <- pair[1]
  #  var2 <- pair[2]
  #  pdf(file.path(fit_name, paste0("joint_density_", fit_name, "_", var1, "_vs_", var2, ".pdf")), width = 5, height = 4)
  #  p_pair <- ggplot(draws_df, aes(x = .data[[var1]], y = .data[[var2]])) + 
  #    geom_density_2d_filled() + 
  #    theme_minimal() + 
  #    labs(x = var1, y = var2, fill = "Density Level")
  #  print(p_pair)
  #  dev.off()
  #}
  
  # 7. Pairs plots
  
  my_diag <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    
    q <- quantile(x, c(0.025, 0.975))
    mu <- mean(x)
    
    lines_df <- data.frame(
      xintercept = c(mu, q),
      type = c("Mean", "95% CI", "95% CI")
    )
    
    ggplot(data, mapping) +
      geom_histogram(aes(y = after_stat(density)), bins = 50, fill="lightblue", color="black", alpha=0.5) +
      geom_density(aes(color="Density", linetype="Density"), linewidth=1, key_glyph = "path") +
      geom_vline(data = lines_df, aes(
        xintercept = xintercept,
        color = type,
        linetype = type
      ), linewidth = 0.8 ) +
      scale_color_manual(name = NULL, 
                         values = c("Density" = "blue", "Mean" = "red", "95% CI" = "green")) +
      scale_linetype_manual(name = NULL, 
                            values = c("Density" = "solid", "Mean" = "solid", "95% CI" = "dashed")) +
      theme_minimal() +
      theme(
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
        legend.text = element_text(size = 8),
        legend.key.width = unit(1, "cm")
      )
  }
  
  p <- ggpairs(
    draws_df[, param_names],
    upper = list(continuous = wrap("density", alpha = 0.7)),
    lower = list(continuous = wrap("cor", stars = FALSE, size = 5)),
    diag = list(continuous = my_diag)
  )
  
  ggsave(
    file.path(fit_name, paste0("pairs_", fit_name, ".pdf")),
    p,
    width = 10,
    height = 10
  )
  print(paste0("Plots saved in directory: ", fit_name))
  print(paste0("End of MCMC analysis"))
  return(summary_stats)
}

stan_data <- list(
  N = nrow(annual),
  t = annual$year - min(annual$year),
  y = annual$annual_mean
)

linear_model <- cmdstan_model("linear_regression.stan")
linear_fit <- linear_model$sample(data=stan_data,
                                  chains=4, 
                                  parallel_chains=4, 
                                  iter_warmup=2000, 
                                  iter_sampling=10000, 
                                  seed=42,
                                  refresh = 0
)
linear_fit$cmdstan_diagnose()

# Posterior Analysis
linear_summary <- analyze_MCMC(
  fit = linear_fit,
  model = linear_model,
  stan_data = stan_data,
  param_names = c("b0", "b1", "sigma"),
  fit_name = "linear"
)



# ================================================================ #
#                      QUADRATIC MODEL FITTING                     #
# ================================================================ #

quadratic_model <- cmdstan_model("quadratic_regression.stan")
quadratic_fit <- quadratic_model$sample(data=stan_data,
                                        chains=4, 
                                        parallel_chains=4, 
                                        iter_warmup=2000, 
                                        iter_sampling=10000,
                                        #thin = 10,
                                        seed=42,
                                        refresh = 0
)
quadratic_fit$cmdstan_diagnose()

# Posterior Analysis
quadratic_summary <- analyze_MCMC(
  fit = quadratic_fit,
  model = quadratic_model,
  stan_data = stan_data,
  param_names = c("b0", "b1", "b2", "sigma"),
  fit_name = "quadratic"
)


# ================================================================ #
#                      MODEL COMPARISON                            #
# ================================================================ #

# comparison using LOO
library(loo)
log_lik_linear <- linear_fit$draws("log_lik")
log_lik_quadratic <- quadratic_fit$draws("log_lik")

# compute LOO for both models
loo_linear <- loo(log_lik_linear)
loo_quadratic <- loo(log_lik_quadratic)
comparison <- loo_compare(loo_linear, loo_quadratic)
print(comparison)


# comparison using Bayes factor ???


# ================================================================ #
#      REGRESSION CURVE AND POSTERIOR PREDICTIVE INTERVAL          #
# ================================================================ #

plot_posterior_predictive <- function(fit, model_name, stan_data, annual_df) {
  # model_name must be either "linear" or "quadratic"
  
  # 1. Dynamically select parameters based on the model
  is_quad <- model_name == "quadratic"
  params <- if (model_name == "quadratic") c("b0", "b1", "b2", "sigma") else c("b0", "b1", "sigma")
  
  draws <- as_draws_df(fit$draws(variables = params))
  
  # 2. Define the time axis
  start_year <- min(annual_df$year)
  t_seq <- seq(min(stan_data$t), max(stan_data$t) + 10, length.out = 1000) 
  
  # 3. Compute the mean curve for each MCMC sample
  # We get a (n_samples x n_time_points) matrix of mu values
  if (is_quad) {
    mu_matrix <- sapply(t_seq, function(t) draws$b0 + draws$b1 * t + draws$b2 * (t^2))
  } else {
    mu_matrix <- sapply(t_seq, function(t) draws$b0 + draws$b1 * t)
  }
  
  # 4. Compute mean curve credible interval
  # We calculate the mean over the columns (i.e., over the samples) for each time point
  curve_mean <- apply(mu_matrix, 2, mean)
  curve_lwr  <- apply(mu_matrix, 2, quantile, probs = 0.025)
  curve_upr  <- apply(mu_matrix, 2, quantile, probs = 0.975)
  
  # 5. Compute posterior predictive interval
  y_pred_matrix <- mu_matrix + matrix(rnorm(length(mu_matrix), mean = 0, sd = draws$sigma), 
                                      nrow = nrow(mu_matrix), ncol = ncol(mu_matrix))
  pred_lwr <- apply(y_pred_matrix, 2, quantile, probs = 0.025)
  pred_upr <- apply(y_pred_matrix, 2, quantile, probs = 0.975)
  
  # 6. Build the consolidated dataframe for plotting
  fit_df <- data.frame(
    year      = t_seq + start_year,
    mu        = curve_mean,
    curve_lwr = curve_lwr,
    curve_upr = curve_upr,
    pred_lwr  = pred_lwr,
    pred_upr  = pred_upr
  )
  
  # 6. Generate and save the plot
  pdf_name <- paste0(model_name, "_fit_with_intervals.pdf")
  pdf(pdf_name, width = 8, height = 5)
  
  p <- ggplot() +
    geom_ribbon(data = fit_df, aes(x = year, ymin = pred_lwr, ymax = pred_upr, fill = "95% Prediction Interval"), alpha = 0.2) +
    geom_ribbon(data = fit_df, aes(x = year, ymin = curve_lwr, ymax = curve_upr, fill = "95% Credible Interval"), alpha = 0.5) +
    geom_line(data = fit_df, aes(x = year, y = mu, color = "Mean Trend"), linewidth = 0.5) +
    geom_point(data = annual_df, aes(x = year, y = annual_mean), color = "black", size = 1, alpha = 0.7) +
    geom_vline(
      xintercept = max(annual_df$year),
      linetype = "dashed",
      color = "black",
      alpha = 0.6
    )+
    scale_fill_manual(name = NULL, values = c("95% Credible Interval" = "red", "95% Prediction Interval" = "grey50")) +
    scale_color_manual(name = NULL, values = c("Mean Trend" = "darkred")) +
    scale_x_continuous(
      breaks = pretty(fit_df$year, 10),
      expand = expansion(mult = c(0.01, 0.05))
    ) +
    scale_y_continuous(
      breaks = pretty(fit_df$mu, 8),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    labs(
      x = "Year",
      y = expression(paste("CO" ["2"], " (ppm)")),
      title = paste("Posterior fit -", model_name, "model")
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      text = element_text(size = 12),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
  print(p)
  dev.off()
  
  cat(sprintf("Plot successfully saved to: %s\n", pdf_name))
  
  future_years <- seq(max(annual_df$year) + 1,
                      max(annual_df$year) + 10,
                      by = 1)
  
  future_df <- subset(fit_df, year %in% future_years)
  
  interval_table <- data.frame(
    year = future_df$year,
    mean = future_df$mu,
    ci_lower = future_df$curve_lwr,
    ci_upper = future_df$curve_upr
  )
  
  print(interval_table)
}


# Linear model
plot_posterior_predictive(fit = linear_fit, 
                          model_name = "linear", 
                          stan_data = stan_data, 
                          annual_df = annual)

# Quadratic model
plot_posterior_predictive(fit = quadratic_fit, 
                          model_name = "quadratic", 
                          stan_data = stan_data, 
                          annual_df = annual)
