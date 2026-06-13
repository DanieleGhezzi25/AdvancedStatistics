data {
    int<lower=1> N;
    vector[N] t;
    vector[N] y; 
}

parameters {
    real<lower=0, upper=1000> b0; // they have to be restricted to uniform prior's values      
    real<lower=-1.56, upper=1.56> alpha;  // not pi/2 exactly otherwise it diverges
    real<lower=-10, upper=2> logsigma;
}

transformed parameters {
    vector[N] mu;
    real<lower=0> sigma;
    real b1; 
    
    b1 = tan(alpha);
    mu = b0 + b1 * t;
    sigma = exp(logsigma);
}

model {
    // Priors
    b0 ~ uniform(0, 1000);
    alpha ~ uniform(-1.56, 1.56);
    logsigma ~ uniform(-10, 2); // jeffrey's prior 1/sigma, with sigma bounded between [e^-10,e^2]
    
    // Likelihood
    y ~ normal(mu, sigma);
}

generated quantities {
    array[N] real y_rep;
    vector[N] log_lik;    
    for (i in 1:N) {
        y_rep[i] = normal_rng(mu[i], sigma); 
        log_lik[i] = normal_lpdf(y[i] | mu[i], sigma); 
    }
}
