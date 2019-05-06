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
fish = read.csv('grasscarplengths.csv')  
  
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

# Simple VBGF in stan -----
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

  
# Hydrilla VBGF in stan -----
# Package the data for stan
  vb_data = list(
    length = fish$Length,
    age = fish$Age,
    nFish = nrow(fish),
    group = as.numeric(as.factor(fish$yearc)),
    ngroups = length(unique(fish$yearc)),
    hydrilla = as.vector(scale(fish$ha))
  )

# Inits for stan
  inits <- function(){
    list(
      b0_linf = rnorm(1, 1000, 50),
      K = rnorm(1, 0, 1),
      t0 = rnorm(1, 0, 1)
    )
  }
  
# Fit the model with stan  
  fit <- stan(file = 'vonbert_hydrilla.stan',
              data = vb_data,
              chains = 3,
              iter = 1000,
              init = inits,
              control = list(adapt_delta = .8)
              )  
  
# Print model summary
  print(fit, digits=3)

# Save fitted object to file
  save(fit, file='hydrillamod-result.rda')
  
  
  
