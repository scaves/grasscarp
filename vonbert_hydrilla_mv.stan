data {
  int<lower = 0>           nobs; // number of fish
  vector<lower = 0>[nobs]  age; // Age of fish
  real<lower = 0>          length[nobs]; // Length of fish
  vector[nobs]             hydrilla; // standarized hydrilla ha
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
  real<lower = 0> sigmaLength;
  real  b0_linf;
  real  bh_linf;
  real  b0_k;
  real  bh_k;
  real  b0_t0;  
}
transformed parameters{
  vector[3] mu_beta_cor = // correlated site-level variation, without mu_gamma
    diag_pre_multiply(tau, L_Omega) * mu_beta_raw;
    // This is part of the MVN non-central parameterization
    // The mean vector (mu_gamma) still needs to be added, but that is done below
  real  Linf = //theoretical maximum length
    exp(mu_beta_cor[1] + b0_linf);
  real  K = // growth coefficient
    exp(mu_beta_cor[2]  + b0_k);
  real  t0 = // hypothetical age at which fish's size = 0
    exp(mu_beta_cor[3] + b0_t0) - 10;
}
model {
  vector[nobs] vonBplaceholder;
  real linf;
  real k;
  
  // Model for expected length at age  
  for(fish in 1:nobs){
  
    linf = Linf + bh_linf * hydrilla[fish];
    k = K + bh_k * hydrilla[fish];
  
    vonBplaceholder[fish] = linf * (1 - exp( -k * ( age[fish] - t0 )));
    
  }
  
  L_Omega ~ lkj_corr_cholesky(hp_omega);
  to_vector(mu_beta_raw) ~ normal(0,1);
  tau ~ exponential(1/hp_tau);
  sigmaLength ~ exponential(1/hp_sigma);

  // estimation
  length ~ normal(vonBplaceholder, sigmaLength);
  
}
