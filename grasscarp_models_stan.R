# Front-end needs ----
# . Load necessary packages ----
  library(Rcpp)
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
fish = read.csv('data/grasscarplengths.csv', stringsAsFactors = F)  
  
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
  hydrilla = read.csv('data/gcStockingAndHydrilla.csv')

# Change year to match names in fish
  colnames(hydrilla)[1] = 'Year'

# Merge hydrilla data with fish data
  fish = merge(x = fish, y = hydrilla, by = 'Year')

# Annual VBGF -----
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
  fish = fish[fish$agec == fish$Age, ]  
  
# Package the data for stan
  vb_data = list(
    length = fish$Length,
    age = fish$Age,
    nobs = nrow(fish),
    fish = sort(unique(as.numeric(as.factor(fish$fishID)))),
    nfish = length(unique(as.numeric(as.factor(fish$fishID)))),
    fishID = as.numeric(as.factor(fish$fishID)),
    hydrilla = as.vector(scale(fish$ha)),
    hp_tau = 1.5,
    hp_sigma = 10,
    hp_omega = 1, # default was 2
    p_linf = 0,
    p_linf_sd = 1,
    p_k = 0,
    p_k_sd = 1
  )

# Parameters to save ---
  params = c('b0_linf', 'bh_linf', 'b0_k', 'bh_k', 'b0_t0', 'mu_beta_cor')
  
# Inits for stan
  inits <- function(){
    list(
      b0_linf = rnorm(1,0,1),
      b0_k = rnorm(1,0,1),
      b0_t0 = rnorm(1,0,1)
    )
  }
  
# Fit the model with stan  
  fit <- stan(file = 'models/vonbert_hydrilla_mv.stan',
              data = vb_data,
              pars = params,
              chains = 4,
              iter = 300,
              warmup = 250,
              init = inits,
              control = list(adapt_delta = .90,
                             max_treedepth=15)
              )  
  
# Print model summary
  print(fit, digits=3)

# Extract parameters
  pars <- extract(fit)
  
# Look at pairs plot
  #pairs(fit)
  
# Sanity check
  hist(exp(pars$b0_k))
  hist(exp(pars$b0_linf))
  
  hist(pars$b0_linf)  
  hist(pars$bh_linf)  
  hist(pars$b0_k)
  hist(pars$bh_k)
  hist(pars$b0_t0)
  
# Save the model fit 
save(fit, file='results/vonbert_hydrilla_mv.rda')  
  

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
  fit <- stan(file = 'models/vonbert_individual_mv.stan',
              data = vb_data,
              pars = params,
              chains = 3,
              iter = 6000,
              warmup = 5000,
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
#  fish=fish[fish$agec >= 15, ]
  
# Package the data for stan
  vb_data = list(
    length = fish$Length,
    age = fish$Age,
    obs = as.numeric(row.names(fish)),
    nobs = nrow(fish),
    fish = sort(unique(as.numeric(as.factor(fish$fishID)))),
    nfish = length(unique(as.numeric(as.factor(fish$fishID)))),
    fishID = as.numeric(as.factor(fish$fishID)),
    hydrilla = as.vector(scale(fish$ha)),
    hp_tau = 1.5,
    hp_sigma = 10,
    hp_omega = 1, # default was 2
    p_linf = 7,
    p_linf_sd = .5,
    p_k = -1,
    p_k_sd = .5
  )

# Parameters to save ---
  params = c('b0_linf', 'bh_linf', 'b0_k', 'bh_k', 'b0_t0')
  
# Inits for stan
  inits <- function(){
    list(
      linf = 7,
      k = -1.4,
      b0_t0 = 0
    )
  }
  
# Fit the model with stan  
  fit <- stan(file = 'vonbert_hydrilla_mv_ind.stan',
              data = vb_data,
              pars = params,
              chains = 3,
              iter = 1000,
              warmup = 800,
              #init = inits,
              control = list(adapt_delta = .80,
                             max_treedepth=10)
              )  
  
# Print model summary
  print(fit, digits=3)

# Extract parameters
  pars <- extract(fit)
  
# Look at pairs plot
  pairs(fit)
  
# Sanity check
  hist(exp(pars$b0_k))
  
  
  
# Save the model fit 
save(fit, file='vonbert_hydrilla_mv_ind.rda')  




# For Jess ----

par(mfrow=c(2, 1), oma=c(4,4,0,1), mar=c(1,1,1,1))
hist(fish2$Length[fish2$Year<=2010],
     main = '', col='gray40',
     yaxt = 'n', ylim = c(0, 60),
     xaxt = 'n', xlim = c(0, 1400),
     breaks = 20
     )
text(x=100, y=50, '2006 - 2010', adj = 0)
axis(1, labels=F, pos=0)
axis(2, las=2, pos=0)
hist(fish2$Length[fish2$Year==2017],
     main = '', col='gray40',
     yaxt = 'n', ylim = c(0, 60),
     xaxt = 'n', xlim = c(0, 1400)
     )
text(x=100, y=50, '2017', adj = 0)
axis(1, labels=TRUE, pos=0)
axis(2, las=2, pos=0)
mtext(text = 'Total length (mm)', side = 1, line=3)
mtext(text = 'Frequency', side = 2, line=3, adj=1.5)
