data {
  int<lower = 0>           nobs; // Number of observations
  int<lower = 0>           nfish; // Number of fish
  int<lower = 0>           fishID[nobs] ; // Fish ID
  row_vector<lower = 0>[nobs] age; // Age of fish
  row_vector<lower = 0>[nobs] length; // Length of fish
  row_vector[nobs] hydrilla; // Standardized hydrilla ha
 // Hyperparameters below:
  real hp_tau;  // default 2.5
  real hp_sigma;  // unsure what default is
  real hp_omega; // default 2
  // prior for mu_gamma
  real p_mu_gamma;
  real p_mu_gammaSD;
}
parameters {
  cholesky_factor_corr[3] L_Omega; // prior correlation, Cholesky scale
  vector[3] mu_beta_raw; // This will be transformed into mu_beta
  vector<lower=0>[3] tau; // prior scale
  real<lower = 0> sigma; // observation error
  real  b0_linf;
  real  bh_linf;
  real  b0_k;
  real  bh_k;
  real  b0_t0;  
}
transformed parameters {
  row_vector[3] mu_gamma[nobs];
  vector[3] mu_beta_cor;
  row_vector[3] theta[nobs];
  mu_gamma[1] = b0_linf + bh_linf*hydrilla;
  mu_gamma[2] = b0_k + bh_k*hydrilla;
  mu_gamma[3] = b0_t0;
  mu_beta_cor = // correlated site-level variation, without mu_gamma
    diag_pre_multiply(tau, L_Omega) * mu_beta_raw;
  theta[1] = //theoretical maximum length
    exp(mu_beta_cor[1] + mu_gamma[1]);
  theta[2] = // growth coefficient
    exp(mu_beta_cor[2]  + mu_gamma[2]);
  theta[3] = // hypothetical age at which fish's size = 0
    exp(mu_beta_cor[3] + mu_gamma[3]) - 10.0;
}
model{
  row_vector[nobs] y;
  
  b0_linf ~ normal(0, 1);
  b0_k ~ normal(0, 1);
  b0_t0 ~ normal(0, 1);
  bh_k ~ normal(0, 1);
  bh_linf ~ normal(0, 1);
  
  L_Omega ~ lkj_corr_cholesky(hp_omega);
  mu_beta_raw ~ normal(0,1);
  tau ~ exponential(1/hp_tau);
  sigma ~ exponential(1/hp_sigma);
  
  y = theta[1] .* (1 - exp( -theta[2] .* ( age - theta[3])));
  
  length ~ normal(y, sigma);
}
