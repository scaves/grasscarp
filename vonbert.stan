data {
  int<lower = 0>           nFish; // number of fish
  int<lower = 0>           ngroups; // number of sites
  int<lower = 0>           group[nFish] ; // Dummy variable to ID pool
  vector<lower = 0>[nFish] age; // Age of fish
  real<lower = 0>          length[ nFish]; // Length of fish
}
parameters {
  real<lower = 0> sigmaLength;
  vector[ngroups]  Linf;
  vector[ngroups]  t0;
  vector[ngroups]  K;
}
model {
  vector[nFish] vonBplaceholder;
  for(fish in 1:nFish){
    vonBplaceholder[fish] = Linf[ group[fish] ] * (1 - exp( -K[group[fish]] * ( age[fish] - t0[group[fish]] )));
  }
  // priors 
  //Linf ~ normal(1200, 50);
  //K    ~ normal(0, 1);
  //t0   ~ normal(0, 1);

  // estimation
  length ~ normal(  vonBplaceholder, sigmaLength);
}
