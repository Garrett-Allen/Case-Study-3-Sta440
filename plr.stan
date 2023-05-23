data {
  int<lower = 1> n; // number of spatial locations
  int<lower = 1> T; // number of longitudinal observations
  int<lower = 1> N; // total number of observations
  vector[N] y; // outcome observations
  vector[N] x; // time measurements
  int<lower = 1> s[N]; // location indeces
}
parameters {
  vector[n] beta0;
  vector[n] beta1;
  vector<lower = 0>[n] sigma;
}
model {
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = beta0[s[i]] + beta1[s[i]] * x[i];
  }
  for (i in 1:n) {
    sigma[i] ~ student_t(3, 0, 1);
  }
  for (i in 1:N) {
    y[i] ~ normal(mu[i], sigma[s[i]]);
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = beta0[s[i]] + beta1[s[i]] * x[i];
  }
  for (i in 1:N) log_lik[i] = normal_lpdf(y[i] | mu[i], sigma[s[i]]);
}
