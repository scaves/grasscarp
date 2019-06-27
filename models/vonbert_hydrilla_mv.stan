data {
  int<lower = 0>           nobs;     // Number of observations
  vector<lower = 0>[nobs]  age;      // Age of fish
  vector<lower = 0>[nobs]  length;   // Length of fish
  vector[nobs]             hydrilla; // Standardized hydrilla ha
  // Priors passed as data
  real p_tau;  
  real p_sigma; 
  real p_omega;
  real p_linf;
  real p_linf_sd;
  real p_k;
  real p_k_sd;  
}
parameters {
  cholesky_factor_corr[3] L_Omega; // Cholesky factor
  vector[3] z;            // Uncorrelated parameter values
  vector<lower=0>[3] tau; // Scale vector
  real<lower = 0> sigma;  // Observation error
  real  b0_linf;          // Intercept for Linf
  real  bh_linf;          // Slope for HA on Linf
  real  b0_k;             // Intercept for K
  real  bh_k;             // Slope for HA on K
  real  b0_t0;            // Mean t0 (shared)
}
transformed parameters {
  vector[3] Gamma;    // Correlated offsets
  vector[nobs] Linf;  // Asymptotic length
  vector[nobs] K;     // Brody growth coeff
  real t0;            // Age at length = 0
  
  Gamma = diag_pre_multiply(tau, L_Omega) * z;
  Linf = exp(Gamma[1] + (b0_linf + bh_linf*hydrilla));
  K = exp(Gamma[2] + (b0_k + bh_k*hydrilla));
  t0 = exp(Gamma[3] + (b0_t0)) - 10;          
}
model{
  vector[nobs] y;
  vector[nobs] vt0 = rep_vector(t0, nobs);
  to_vector(z) ~ normal(0,1);

  y = Linf .* (1-exp(-K .* (age-vt0)));
  
  b0_linf ~ normal(p_linf, p_linf_sd);
  b0_k ~ normal(p_k, p_k_sd);
  b0_t0 ~ normal(0, 1);
  bh_k ~ normal(0, 1);
  bh_linf ~ normal(0, 1);
  
  L_Omega ~ lkj_corr_cholesky(p_omega);
  tau ~ exponential(1/p_tau);
  sigma ~ exponential(1/p_sigma);
  
  length ~ normal(y, sigma);
}
