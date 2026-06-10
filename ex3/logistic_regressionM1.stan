data {
    int<lower=1> N;
    vector[N] x;
    array[N] int<lower=0, upper=1> y;
}
parameters {
    real beta0;
    real beta1;
}
model {
    beta0 ~ normal(0, 3);
    beta1 ~ normal(0, 3);
    y ~ bernoulli_logit(beta0 + beta1 * x);
}
generated quantities {
    array[N] int y_rep;
    vector[N] log_lik;
    for (i in 1:N) {
        y_rep[i] = bernoulli_logit_rng(beta0 + beta1 * x[i]); // samples
        log_lik[i] = bernoulli_logit_lpmf(y[i] | beta0 + beta1 * x[i]); // log pmf
    }
}
