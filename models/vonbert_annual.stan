data {
  int<lower = 0>                  nFish; // number of fish
  int<lower = 0>                  ngroups; // number of sites
  int<lower = 0>                  group[nFish]; // Dummy variable to ID pool
  row_vector<lower = 0>[nFish]    age; // Age of fish
  row_vector<lower = 0>[nFish]    length; // Length of fish
  real                            hp_omega;
  real                            hp_tau;
  real                            hp_sigma;
  real                            p_mu_gamma;
  real                            p_mu_gamma_sd;
  
}
parameters {
  cholesky_factor_corr[3] L_Omega; // prior correlation, Cholesky scale
  matrix[3, ngroups]      mu_beta_raw; // This will be transformed into mu_beta
  vector<lower=0>[3]      tau; // prior scale
  real<lower = 0>         sigma_length; // observation error
  vector[3]               mu_gamma; // group coeffs
}
transformed parameters {
  matrix[3, ngroups]   mu_beta_cor;
  row_vector[ngroups]  Linf;
  row_vector[ngroups]  K;
  row_vector[ngroups]  t0;
  
  mu_beta_cor = diag_pre_multiply(tau, L_Omega) * mu_beta_raw;
  Linf = exp(mu_beta_cor[1] + mu_gamma[1]);
  K = exp(mu_beta_cor[2]  + mu_gamma[2]);
  t0 = exp(mu_beta_cor[3] + mu_gamma[3]) - 10.0;
}
model {
  row_vector[nFish] vb;
  
  vb = Linf[group] .* (1 - exp( -K[group] .* (age - t0[group])));
  
  // estimation
  length ~ normal(vb, sigma_length);
  
  L_Omega ~ lkj_corr_cholesky(hp_omega);
  to_vector(mu_beta_raw) ~ normal(0,1);
  tau ~ exponential(1/hp_tau);
  sigma_length ~ exponential(1/hp_sigma);
  mu_gamma ~ normal( p_mu_gamma, p_mu_gamma_sd);
  
}
