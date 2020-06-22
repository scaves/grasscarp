data {
  int<lower = 0>              nobs;        // Number of observations
  int<lower = 0>              ngroups;     // number of groups
  int<lower = 0>              group[nobs]; // Group ID
  row_vector<lower = 0>[nobs] age;         // Age of fish
  row_vector<lower = 0>[nobs] length;      // Length of fish
  // Priors passed as data
  real p_linf;
  real p_linfsd;
  real p_k;
  real p_ksd;
  real p_t0;
  real p_t0sd;
  real p_sigma;
}
parameters {
  real<lower = 0>             sigma;   // Observation error
  row_vector[ngroups]         linf;    // Global VBGF params
  row_vector[ngroups]         k;       // Global VBGF params
  row_vector[ngroups]         lt0;     // Global VBGF params
  real         linf_mu;    // Global VBGF params
  real         k_mu;       // Global VBGF params
  real         t0_mu;     // Global VBGF params  
}
transformed parameters {
  row_vector[ngroups]  Linf;      // Asymptotic length
  row_vector[ngroups]  K;         // Brody growth coefficient
  row_vector[ngroups]  t0;        // Age at length 0
  
  Linf = exp(linf);
  K = exp(k);
  t0 = exp(lt0) - 10.0;
}
model {
  row_vector[nobs] y;
  
  // Length as expectation of VBGF
  y = Linf[group] .* (1 - exp( -K[group] .* (age - t0[group])));
  
  linf_mu ~ normal(p_linf, p_linfsd);
  k_mu ~ normal(p_k, p_ksd);
  t0_mu ~ normal(p_t0, p_t0sd);
  
  for(i in 1:ngroups){
    linf ~ normal(linf_mu, 1);
    k ~ normal(k_mu, 1);
    lt0 ~ normal(t0_mu, 1);    
  }
  
  //Process error
  sigma ~ exponential(1/p_sigma);
  
  // Likelihood
  length ~ normal(y, sigma);

}
