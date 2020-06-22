data {
  int<lower = 0> nobs;              // Number of observations
  int<lower = 0>  obs[nobs];        // Observation
  row_vector<lower = 0>[nobs] age;  // Age of fish
  vector<lower = 0>[nobs] length;   // Length of fish
  row_vector[nobs] x;               // Standardized abund ha
  int<lower = 0> ngroup;            // Number of cohorts
  int<lower = 0> group[nobs];       // Cohort of fish
  
  // Priors passed as data
  real p_tau;                       // Scale (1)
  real p_sigma;                     // Observation error (10)
  real p_omega;                     // Cholesky prior for LKJ (4)
  real p_linf;                      // Prior on mean b0_linf (0)
  real p_linf_sd;                   // Prior on sd b0_linf (1)
  real p_b_linf_sd;                 // Prior on sd b_linf (1)
  real p_k;                         // Prior on mean b0_k (0)
  real p_k_sd;                      // Prior on sd b_k (1)
  real p_b_k_sd;
  
}
parameters {
  cholesky_factor_corr[3] L_Omega;  // Cholesky factor
  matrix[3, ngroup] z;              // Uncorrelated parameter values
  vector<lower=0>[3] tau;           // Scale vector
  real<lower = 0> sigma;            // Observation error
  real  b0_linf;                    // Intercept for Linf
  real  ba_linf;                    // Slope for X on Linf
  real  b0_k;                       // Intercept for K
  real  ba_k;                       // Slope for X on K
  real  b0_t0;                      // Mean t0 (shared)
  
}
transformed parameters {
  matrix[3, ngroup] Gamma;          // Correlated offsets
  row_vector[nobs] Linf;            // Asymptotic length
  row_vector[nobs] K;               // Brody growth coeff
  row_vector[nobs] t0;              // Age at length = 0
  
  // Correlated offsets to VBGF Parameters
  Gamma = diag_pre_multiply(tau, L_Omega) * z;
  
  // VBGF parameters wtih cohort-specific, correlated offsets
  // and covariate effects on Linf and K
  Linf = exp(Gamma[1, group] + (b0_linf + ba_linf * x[obs]));
  K = exp(Gamma[2, group] + (b0_k + ba_k*x[obs]));
  t0 = exp(Gamma[3, group] + (b0_t0)) - 10;   
  
}
model{
  row_vector[nobs] y;
  to_vector(z) ~ normal(0,1);

  // Log-scale priors on VBGF parameters
  b0_linf ~ normal(p_linf, p_linf_sd);
  b0_k ~ normal(p_k, p_k_sd);
  b0_t0 ~ normal(0, 1);
  ba_k ~ normal(0, p_b_k_sd);
  ba_linf ~ normal(0, p_b_linf_sd);
  
  // Priors on VBGF parameter correlation
  L_Omega ~ lkj_corr_cholesky(p_omega);
  tau ~ exponential(1/p_tau);
  
  // Prior for observation error
  sigma ~ exponential(1/p_sigma);
  
  // Expectation of VBGF for each individual
  y = Linf .* (1-exp(-K .* (age-t0)));

  // Likelihood
  length ~ normal(y, sigma);
}
generated quantities{
 vector[nobs] log_lik;
 
 // Likelihood for each HMC iteration
 for(i in 1:nobs) log_lik[i] = normal_lpdf(length[i] | Linf[i] * (1-exp(-K[i]*(age[i]-t0[i]))), sigma);
 
}
