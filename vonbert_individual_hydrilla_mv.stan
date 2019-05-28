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
  matrix[3, nfish] mu_beta_raw; // This will be transformed into mu_beta
  vector<lower=0>[3] tau; // prior scale
  real<lower = 0> sigmaLength; // observation error
  real  b0_linf;
  real  bh_linf;
  real  b0_k;
  real  bh_k;
  real  b0_t0;  
  
}
transformed parameters {
  matrix[3, nfish] mu_beta_cor = // correlated site-level variation, without mu_gamma
    diag_pre_multiply(tau, L_Omega) * mu_beta_raw;
    // This is part of the MVN non-central parameterization
    // The mean vector (mu_gamma) still needs to be added, but that is done below
  row_vector[nfish]  Linf = //theoretical maximum length
    exp(mu_beta_cor[1] + b0_linf);
  row_vector[nfish]  K = // growth coefficient
    exp(mu_beta_cor[2]  + b0_k);
  row_vector[nfish]  t0 = // hypothetical age at which fish's size = 0
    exp(mu_beta_cor[3] + b0_t0) - 10.0;
}
model{
  row_vector[nobs] vonBplaceholder;
  row_vector[nobs] linf;
  vector[nobs] k;
  
  b0_linf ~ normal(0, 1);
  b0_k ~ normal(0, 1);
  b0_t0 ~ normal(0, 1);
  bh_k ~ normal(0, 1);
  bh_linf ~ normal(0, 1);
  
  // Model for expected length at age  
  for(fish in 1:nobs){
  
    linf[fish] = Linf[fishID[fish]] + bh_linf * hydrilla[fish];
    k[fish] = K[fishID[fish]] + bh_k * hydrilla[fish];
  
    vonBplaceholder[fish] = linf[fish]* (1 - exp( -k[fish] * ( age[fish] - t0[fishID[fish]] )));
    
  }
  
  L_Omega ~ lkj_corr_cholesky(hp_omega);
  to_vector(mu_beta_raw) ~ normal(0,1);
  tau ~ exponential(1/hp_tau);
  sigmaLength ~ exponential(1/hp_sigma);
  
  length ~ normal(vonBplaceholder, sigmaLength);
}
