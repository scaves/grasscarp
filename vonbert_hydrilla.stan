data {
  int<lower = 0>           nFish; // number of fish
  int<lower = 0>           ngroups; // number of sites
  int<lower = 0>           group[nFish] ; // Dummy variable to ID pool
  vector<lower = 0>[nFish] age; // Age of fish
  real<lower = 0>          length[nFish]; // Length of fish
  vector[nFish]            hydrilla; // standarized hydrilla ha
  
}
parameters {
  real<lower = 0> sigmaLength;
  //real  linf;
  //real  k;
  real  b0_linf;
  real  bh_linf;
  real  b0_k;
  real  bh_k;
  real  t0;
  //real<lower = 0> s_linf;
  //real<lower = 0> s_k;
  
}
model {
  vector[nFish] vonBplaceholder;
  real linf;
  real k;

  // Priors on linear models of parameters
  b0_linf ~ normal(1100, 100);  
  
  // Model for expected length at age  
  for(fish in 1:nFish){
  
    linf = b0_linf + bh_linf * hydrilla[fish];
    k = b0_k + bh_k * hydrilla[fish];
  
    vonBplaceholder[fish] = linf * (1 - exp( -k * ( age[fish] - t0 )));
    
  }
  
  // estimation
  length ~ normal(vonBplaceholder, sigmaLength);
  
}
