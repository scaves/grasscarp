# Front-end needs ----
# . Load necessary packages ----
  library(rstan)
  library(FSA)
  library(plyr)

# . Set options ----
  options(mc.cores = parallel::detectCores())  
  rstan_options(auto_write = TRUE)
  Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

# . Function definition ----
# ses
  ses <- function(x){
    sd(x)/sqrt(length(x)-1)
  }
  
# Data manipulation ----
# . Fish data ----
fish = read.csv('grasscarplengths.csv', stringsAsFactors = F)  
  
# Drop missing data   
fish <- fish[!is.na(fish$Age) & !is.na(fish$Length), ]

# Change name for year of capture
names(fish)[2] <- 'yearc'

# Get age of capture by aggregation and merge
maxes <- aggregate(Age~fishID, fish, max)
names(maxes)[2] <- 'agec'
fish <- merge(fish, maxes)

# Calculate year for each growth increment
fish$Year <- fish$yearc-(fish$agec-fish$Age)

# Drop highly questionable samples0
fish <- fish[!(fish$yearc==2017 & fish$agec > 21 & fish$Length < 951), ]


# . Hydrilla data -----
# Read data
  hydrilla = read.csv('gcStockingAndHydrilla.csv')

# Change year to match names in fish
  colnames(hydrilla)[1] = 'Year'

# Merge hydrilla data with fish data
  fish = merge(x = hydrilla, y = fish, by = 'Year')

# Simple VBGF -----
# Package the data for stan
  vb_data = list(
    length = fish$Length,
    age = fish$Age,
    nFish = nrow(fish),
    group = as.numeric(as.factor(fish$yearc)),
    ngroups = length(unique(fish$yearc))
  )

# Initial values for parameters
# Get observed max lengths for each 
# age class
  starts <- ddply(fish,
        'yearc',
        summarize,
        linfs = max(Length),
        se=ses(Length)
        )

# Inits for stan
  inits <- function(){
    list(
      Linf = rnorm(length(unique(fish$yearc)), starts$linfs, starts$se),
      b0_k = rnorm(length(unique(fish$yearc)), 0, 1),
      t0 = runif(length(unique(fish$yearc)), -10, 0)
    )
  }
  
# # Fit the model with stan  
#   fit <- stan(file = 'vonbert.stan',
#               data = vb_data,
#               chains = 3,
#               iter = 1000,
#               init = inits
#               )  
# 
# # Print model summary  
#   print(fit, digits=3)
#   
# # Save result to a file
#   save(fit, file='yearmod-result.rda')

  
# Multivariate with hydrilla effect -----
# Package the data for stan
  vb_data = list(
    length = fish$Length,
    age = fish$Age,
    nobs = nrow(fish),
    nfish = length(unique(as.numeric(as.factor(fish$fishID)))),
    fishID = as.numeric(as.factor(fish$fishID)),
    hydrilla = as.vector(scale(fish$ha)),
    hp_tau = 1.5,
    hp_sigma = 10,
    hp_omega = 1, # default was 2
    p_mu_gamma = 0,
    p_mu_gammaSD = 2
  )

# Parameters to save ---
  params = c('mu_gamma', 'beta_linf', 'beta_k')
  
# Inits for stan
  inits <- function(){
    list(
      b0_linf = 1100,
      b0_k = .5,
      b0_t0 = -1
    )
  }
  
# Fit the model with stan  
  fit <- stan(file = 'vonbert_hydrilla_mv.stan',
              data = vb_data,
              #pars = params,
              chains = 3,
              iter = 1000,
              warmup = 900,
              #init = inits,
              control = list(adapt_delta = .8)
              )  
  
# Print model summary
  print(fit, digits=3)

# Extract parameters
  pars <- extract(fit)
  
# Make a histogram of Linf  
  hist(pars$b0_linf) 
  boxplot(pars$Linf[,1:10])
  

# Individual random effect -----
# Package the data for stan
  vb_data = list(
    length = fish$Length,
    age = fish$Age,
    nobs = nrow(fish),
    nfish = length(unique(as.numeric(as.factor(fish$fishID)))),
    fishID = as.numeric(as.factor(fish$fishID)),
    hydrilla = as.vector(scale(fish$ha)),
    hp_tau = 1.5,
    hp_sigma = 10,
    hp_omega = 1, # default was 2
    p_mu_gamma = 0,
    p_mu_gammaSD = 2
  )

# Parameters to save ---
  params = c('mu_gamma', 'Linf', 'K', 't0')
  
# Inits for stan
  inits <- function(){
    list(
      mu_gamma = c(7, 0, 0)
    )
  }
  
# Fit the model with stan  
  fit <- stan(file = 'vonbert_individual_mv.stan',
              data = vb_data,
              pars = params,
              chains = 3,
              iter = 4000,
              warmup = 3500,
              #init = inits,
              control = list(adapt_delta = .8)
              )  
  
# Print model summary
  print(fit, digits=3)

# Extract parameters
  pars <- extract(fit)
  
# Make a histogram of Linf  
  hist(exp(pars$mu_gamma[,1])) 
  boxplot(pars$Linf[,1:100])  
  

# Individual + hydrilla mv STILL WORKING ON THIS -----
# Package the data for stan
  vb_data = list(
    length = fish$Length,
    age = fish$Age,
    nobs = nrow(fish),
    nfish = length(unique(as.numeric(as.factor(fish$fishID)))),
    fishID = as.numeric(as.factor(fish$fishID)),
    hydrilla = as.vector(scale(fish$ha)),
    hp_tau = 1.5,
    hp_sigma = 10,
    hp_omega = 1, # default was 2
    p_mu_gamma = 0,
    p_mu_gammaSD = 2
  )

# Parameters to save ---
  params = c('b0_linf', 'b0_k', 'b0_t0', 'bh_k', 'bh_linf')
  
# Inits for stan
  inits <- function(){
    list(
      b0_linf = 7.1,
      bh_linf = -2
    )
  }
  
# Fit the model with stan  
  fit <- stan(file = 'vonbert_individual_hydrilla_mv.stan',
              data = vb_data,
              pars = params,
              chains = 3,
              iter = 4000,
              warmup = 3500,
              #init = inits,
              control = list(adapt_delta = .8)
              )  
  
# Print model summary
  print(fit, digits=3)

# Extract parameters
  pars <- extract(fit)
  
# Sanity check
  hist(exp(pars$b0_k))
  
# Save the model fit 
save(fit, file='vonbert_individual_hydrilla_mv.rda')  
