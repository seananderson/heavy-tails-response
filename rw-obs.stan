data {
  int<lower=0> N; // rows of data
  vector[N] y; // vector to hold observations
  real<lower=0> sigma_obs;
}
  parameters {
  real<lower=0> sigma_proc;
  vector[N] U; // states
}
  model {
  sigma_proc ~ cauchy(0, 2.5);
  for (i in 2:N) {
    U[i] ~ normal(U[i-1], sigma_proc);
  }
  y ~ normal(U, sigma_obs);
}
