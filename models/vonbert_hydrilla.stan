data {
  int<lower = 0>           nobs; // number of fish
  int<lower = 0>           nfish; // number of sites
  vector<lower = 0>[nobs]  age; // Age of fish
  real<lower = 0>          length[nobs]; // Length of fish
  vector[nobs]             hydrilla; // standarized hydrilla ha
  
}
parameters {
  real<lower = 0> sigmaLength;
  real  b0_linf;
  real  bh_linf;
  real  b0_k;
  real  bh_k;
  real  t0;
  
}
model {
  vector[nobs] vonBplaceholder;
  real linf;
  real k;

  // Priors on linear models of parameters
  b0_linf ~ normal(1100, 100);  
  
  // Model for expected length at age  
  for(fish in 1:nobs){
  
    linf = b0_linf + bh_linf * hydrilla[fish];
    k = b0_k + bh_k * hydrilla[fish];
  
    vonBplaceholder[fish] = linf * (1 - exp( -k * ( age[fish] - t0 )));
    
  }
  
  // estimation
  length ~ normal(vonBplaceholder, sigmaLength);
  
}
