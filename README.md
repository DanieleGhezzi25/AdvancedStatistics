# Advanced Statistics for Physics Analysis
This repository contains the final project and assignments for the **Advanced Statistics for Physics Analysis** course (Physics of Data, University of Padua). 

## Project
This study investigates the atmospheric $CO_2$ concentration records from the Mauna Loa Observatory (1958–2025) to model the dynamics of carbon accumulation using **Bayesian regression**. By comparing linear (constant growth) and quadratic (accelerating growth) empirical models, we provide decisive statistical evidence favoring an accelerating trend in $CO_2$ accumulation. The model projects the underlying atmospheric $CO_2$ trend to reach approximately 453.0 ppm by 2035.

The analysis is conducted within a rigorous Bayesian framework, utilizing the following techniques:
- **Bayesian Inference & Modeling**: Implementation of weakly informative priors and covariate standardization to optimize polynomial model exploration.
- **MCMC Sampling**: Hamiltonian Monte Carlo (HMC) employing the No-U-Turn Sampler (NUTS) via Stan.
- **Model Comparison**: Expected Log Predictive Density (ELPD) evaluation through Leave-One-Out Cross-Validation (LOO-CV) and Pareto Smoothed Importance Sampling (PSIS).
- **Diagnostics & Validation**: Rigorous MCMC convergence checks (Gelman-Rubin $\hat{R}$, ESS) and model validation via Posterior Predictive Checks.
- **Probabilistic Forecasting**: Computation of 95% credible intervals for the mean structural trend and posterior predictive intervals for future observations.

## Assignments
The Assignments directory contains a collection of graded coursework focusing on core statistical techniques, including Bayesian inference, MCMC methods, and Bayesian regression.
