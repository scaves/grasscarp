data {
  int<lower = 0>              nobs;        // Number of observations
  int<lower = 0>              ngroups;     // number of groups
  int<lower = 0>              group[nobs]; // Group ID
  row_vector<lower = 0>[nobs] age;         // Age of fish
  row_vector<lower = 0>[nobs] length;      // Length of fish
  // Priors passed as data
  real hp_omega;
  real hp_tau;
  real hp_sigma;
  real p_b;
  real p_b_sd;
}
parameters {
  cholesky_factor_corr[3] L_Omega; // Cholesky factor
  matrix[3, ngroups]      Z;       // Uncorrelated parameters
  vector<lower=0>[3]      tau;     // Scale
  real<lower = 0>         sigma;   // Observation error
  vector[3]               b0;      // Global VBGF params
}
transformed parameters {
  matrix[3, ngroups]   Gamma;     // Group-specific correlated offsets
  row_vector[ngroups]  Linf;      // Asymptotic length
  row_vector[ngroups]  K;         // Brody growth coefficient
  row_vector[ngroups]  t0;        // Age at length 0
  
  Gamma = diag_pre_multiply(tau, L_Omega) * Z;
  Linf = exp(Gamma[1] + b0[1]);
  K = exp(Gamma[2]  + b0[2]);
  t0 = exp(Gamma[3] + b0[3]) - 10.0;
}
model {
  row_vector[nobs] y;
  
  // Length as expectation of VBGF
  y = Linf[group] .* (1 - exp( -K[group] .* (age - t0[group])));
  
  // Likelihood
  length ~ normal(y, sigma);
  
  L_Omega ~ lkj_corr_cholesky(hp_omega);
  to_vector(Z) ~ normal(0,1);
  tau ~ exponential(1/hp_tau);
  sigma ~ exponential(1/hp_sigma);
  b0 ~ normal(p_b, p_b_sd);
}
