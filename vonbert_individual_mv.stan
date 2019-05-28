data {
  int<lower = 0>           nobs; // number of fish
  int<lower = 0>           nfish; // number of sites
  int<lower = 0>           fishID[nobs] ; // Dummy variable to ID pool
  row_vector<lower = 0>[nobs] age; // Age of fish
  row_vector<lower = 0>[nobs] length; // Length of fish
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
  matrix[3, nfish] mu_beta_raw; // This will be transformed into mu_beta
  vector<lower=0>[3] tau; // prior scale
  real<lower = 0> sigmaLength; // observation error
  vector[3] mu_gamma; // group coeffs
}
transformed parameters {
  matrix[3, nfish] mu_beta_cor = // correlated site-level variation, without mu_gamma
    diag_pre_multiply(tau, L_Omega) * mu_beta_raw;
    // This is part of the MVN non-central parameterization
    // The mean vector (mu_gamma) still needs to be added, but that is done below
  row_vector[nfish]  Linf = //theoretical maximum length
    exp(mu_beta_cor[1] + mu_gamma[1]);
  row_vector[nfish]  K = // growth coefficient
    exp(mu_beta_cor[2]  + mu_gamma[2]);
  row_vector[nfish]  t0 = // hypothetical age at which fish's size = 0
    exp(mu_beta_cor[3] + mu_gamma[3]) - 10.0;
}
model{
  row_vector[nobs] vonBplaceholder =
    Linf[fishID]  .*  (1.0 - exp( - K[fishID] .* ( age - t0[fishID])));
  L_Omega ~ lkj_corr_cholesky(hp_omega);
  to_vector(mu_beta_raw) ~ normal(0,1);
  tau ~ exponential(1/hp_tau);
  sigmaLength ~ exponential(1/hp_sigma);
  mu_gamma ~ normal(p_mu_gamma, p_mu_gammaSD);
  length ~ normal(vonBplaceholder, sigmaLength);
}
