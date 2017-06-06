data {
  int N; // number of samples
  int G; // number of features

  vector [N] y[G]; // gene expression input
  real x[N]; // covariate measurements

  
}

parameters {
  real beta[G];
  real alpha[G];
  // real mu[G];
  real<lower = 0> chi[G]; // ARD parameters on beta
  real<lower = 0> tau[G];
  real c[G];

  real z[N]; // pseudotimes
}

model {
  z ~ normal(0, 1);
  c ~ normal(0, 1);
  alpha ~ normal(0, 1);
  // mu ~ normal(0, 1);
  chi ~ gamma(10, 1);
  tau ~ gamma(2, 2);
  
  for(g in 1:G) {
    beta[g] ~ normal(0, 1 / sqrt(chi[g]));
    for(i in 1:N) {
      y[g,i] ~ normal(alpha[g] * x[i] + (c[g] + beta[g] * x[i]) * z[i], 1 / sqrt(tau[g]));
    }
  }
}

