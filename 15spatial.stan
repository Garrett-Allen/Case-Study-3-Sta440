data {
  int<lower = 1> n; // number of spatial locations
  int<lower = 1> T; // number of longitudinal observations
  int<lower = 1> N; // total number of observations
  vector[N] y; // outcome observations
  vector[N] x; // time measurements
  int<lower = 1> s[N]; // location indeces
  matrix<lower = 0, upper = 1>[n, n] W; // adjacency matrix
}
transformed data {
  vector[n] zeros;
  matrix[n, n] identity;
  matrix<lower = 0>[n, n] D;
  vector[n] W_rowsums;
  for (i in 1:n) {
    W_rowsums[i] = sum(W[i, ]);
  }
  D = diag_matrix(W_rowsums);
  zeros = rep_vector(0, n);
  identity = diag_matrix(rep_vector(1.0, n));
}
parameters {
  real beta0;
  real beta1;
  vector[n] beta0_vec;
  vector[n] beta1_vec;
  vector<lower = 0>[n] sigma;
  real<lower = 0> tau0;
  real<lower = 0, upper = 1> rho0;
  real<lower = 0> tau1;
  real<lower = 0, upper = 1> rho1;
}
transformed parameters {
  cov_matrix[n] precision0 = (1 / (tau0 * tau0)) * (rho0 * (D - W) + (1 - rho0) * identity);
  cov_matrix[n] precision1 = (1 / (tau1 * tau1)) * (rho1 * (D - W) + (1 - rho1) * identity);
}
model {
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = (beta0 + beta0_vec[s[i]]) + (beta1 + beta1_vec[s[i]]) * x[i];
  }
  beta0_vec ~ multi_normal_prec(zeros, precision0);
  beta1_vec ~ multi_normal_prec(zeros, precision1);
  for (i in 1:n) {
    sigma[i] ~ student_t(3, 0, 1);
  }
  tau0 ~ student_t(3, 0, 1);
  tau1 ~ student_t(3, 0, 1);
  for (i in 1:N) {
    y[i] ~ normal(mu[i], sigma[s[i]]);
  }
  rho0 ~ beta(1,5);
  rho1 ~ beta(1,5);
}
generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = (beta0 + beta0_vec[s[i]]) + (beta1 + beta1_vec[s[i]]) * x[i];
  }
  for (i in 1:N) log_lik[i] = normal_lpdf(y[i] | mu[i], sigma[s[i]]);
}
