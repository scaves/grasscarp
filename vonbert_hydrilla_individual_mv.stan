data {
  int<lower = 0>           nobs; // Number of observations
  int<lower = 0>           nfish; // Number of fish
  int<lower = 0>           fish[nfish];  
  int<lower = 0>           fishID[nobs] ; // Fish ID
  vector<lower = 0>[nobs]  age; // Age of fish
  vector<lower = 0>[nobs]  length; // Length of fish
  vector[nobs] hydrilla; // Standardized hydrilla ha
 // Hyperparameters below:
  real hp_tau;  // default 2.5
  real hp_sigma;  // unsure what default is
  real hp_omega; // default 2
  // prior for mu_gamma
  real p_linf;
  real p_linf_sd;
  real p_k;
  real p_k_sd;  
}
parameters {
  cholesky_factor_corr[3] L_Omega; // prior correlation, Cholesky scale
  vector[3] mu_beta_raw; // This will be transformed into mu_beta
  vector<lower=0>[3] tau; // prior scale
  real<lower = 0> sigma; // observation error
  real  linf;
  vector[nfish]  b0_linf;
  real  bh_linf;
  real  k;
  vector[nfish]  b0_k;
  real  bh_k;
  real  b0_t0;  
}
transformed parameters {
  vector[3] mu_beta_cor;
  vector[nobs] Linf;
  vector[nobs] K;
  real t0;
  
  mu_beta_cor = // correlated site-level variation, without mu_gamma
    diag_pre_multiply(tau, L_Omega) * mu_beta_raw;
  Linf = //theoretical maximum length
    exp(mu_beta_cor[1] + (b0_linf[fishID] + bh_linf*hydrilla));
  K = // growth coefficient
    exp(mu_beta_cor[2] + (b0_k[fishID] + bh_k*hydrilla));
  t0 = // hypothetical age at which fish's size = 0
    mu_beta_cor[3] + (b0_t0);
}
model{
  vector[nobs] y;
  vector[nobs] vt0 = rep_vector(t0, nobs);
  to_vector(mu_beta_raw) ~ normal(0,1);

  y = Linf .* (1-exp(-K .* (age-vt0)));
  
  linf ~ normal(p_linf, p_linf_sd);
  k ~ normal(p_k, p_k_sd);
  
  to_vector(b0_linf) ~ normal(linf, .5);
  to_vector(b0_k) ~ normal(k, .5);  
  b0_t0 ~ normal(0, 1);
  bh_k ~ normal(0, 1);
  bh_linf ~ normal(0, 1);
  
  L_Omega ~ lkj_corr_cholesky(hp_omega);
  tau ~ exponential(1/hp_tau);
  sigma ~ exponential(1/hp_sigma);
  
  length ~ normal(y, sigma);
}
