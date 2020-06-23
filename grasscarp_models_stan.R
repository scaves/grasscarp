# Front-end needs ----
# . Load necessary packages ----
  library(Rcpp)
  library(rstan)
  library(plyr)
  library(loo)

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
fish$cohort <- fish$yearc - fish$agec

# Drop highly questionable samples
fish <- fish[!(fish$yearc==2017 & fish$agec > 21 & fish$Length < 951), ]


# . Hydrilla data -----
# Read data
  hydrilla = read.csv('data/gcStockingAndHydrilla.csv')

# Change year to match names in fish
  colnames(hydrilla)[1] = 'Year'

# Merge hydrilla data with fish data
  fish = merge(x = fish, y = hydrilla, by = 'Year')
  
# Raw data visualization ----
# . Table 2 ----
# Get counts of fish collected each year for
  fishies = fish[fish$agec==fish$Age,]
  totals = ddply(fishies, 'Year', summarize, nums = length(fishID))
  
# . Catch-curve analysis fail ----
  ccs = ddply(fishies, c('Year', 'Age'), summarize, nums = length(fishID))
  ccs17 = ccs[ccs$Year==2017,] 
  summary(lm(nums~Age, data=ccs17[ccs17$Age>=9, ]))
  
# # Annual VBGF -----
#   fish = fish[fish$agec == fish$Age, ]  
#   
# # Package the data for stan
#   vb_data = list(
#     length = fish$Length,
#     age = fish$Age,
#     nobs = nrow(fish),
#     group = as.numeric(as.factor(fish$yearc)),
#     ngroups = length(unique(fish$yearc)),
#     hp_tau = 2.5,
#     hp_sigma = 10,
#     hp_omega = 2,
#     p_b = 0,
#     p_b_sd = 1
#   )
#   
# # Initial values
#   inits = function(){
#     list(
#       b0 = c(7,-2,0)
#     )
#   }
# 
# # Fit the model with stan
#   # fit <- stan(file = 'models/vonbert_annual.stan',
#   #             data = vb_data,
#   #             chains = 3,
#   #             iter = 5000,
#   #             warmup = 4500,
#   #             control = list(
#   #               adapt_delta = 0.999,
#   #               max_treedepth = 15
#   #             )
#   #             )
# 
# # Print model summary
#   print(fit, digits=3)
# 
# # Save result to a file
#   save(fit, file='results/vonbert_annual.rda')
# 
# 
# # Multivariate with hydrilla effect -----
#   fish = fish[fish$agec == fish$Age, ]  
# 
# # Package the data for stan
#   vb_data = list(
#     length = fish$Length,
#     age = fish$Age,
#     nobs = nrow(fish),
#     fish = sort(unique(as.numeric(as.factor(fish$fishID)))),
#     nfish = length(unique(as.numeric(as.factor(fish$fishID)))),
#     fishID = as.numeric(as.factor(fish$fishID)),
#     hydrilla = as.vector(scale(fish$ha)),
#     p_tau = 2.5,
#     p_sigma = 10,
#     p_omega = 2,
#     p_linf = 0,
#     p_linf_sd = 1,
#     p_k = 0,
#     p_k_sd = 1
#   )
# 
# # Fit the model with stan  
#   fit <- stan(file = 'models/vonbert_hydrilla_mv.stan',
#               data = vb_data,
#               chains = 3,
#               iter = 5000,
#               warmup = 4500,
#               control = list(adapt_delta = .80,
#                              max_treedepth=10)
#               )  
#   
# # Print model summary
#   print(fit, digits=3)
#   
# # Save the model fit 
# save(fit, file='results/vonbert_hydrilla_mv.rda')  
#   
# # Cohort VBGF -----
#   #fish = fish[fish$agec == fish$Age, ]  
#   
# # Package the data for stan
#   vb_data = list(
#     length = fish$Length,
#     age = fish$Age,
#     nobs = nrow(fish),
#     group = as.numeric(as.factor(fish$cohort)),
#     ngroups = length(unique(fish$cohort)),
#     p_linf = 0,
#     p_linfsd = 1,
#     p_k = 0,
#     p_ksd = 1,  
#     p_t0 = 0,
#     p_t0sd = 1,      
#     p_sigma = 10
#   )
# 
# # Fit the model with stan
#   fit <- stan(file = 'models/vonbert_cohort.stan',
#               data = vb_data,
#               chains = 3,
#               iter = 3000,
#               #warmup = 4500,
#               control = list(
#                 adapt_delta = 0.80,
#                 max_treedepth = 10
#               )
#               )
# 
# # Print model summary
#   print(fit, digits=3)
# 
# # Save result to a file
#   save(fit, file='results/vonbert_cohort.rda')
# 
# # Extract parameters
#   pars = extract(fit)
#   boxplot(pars$Linf, outline=FALSE)  
#   boxplot(pars$K, outline=FALSE)  
#   
#   
# Cohort x covariate VBGF -----
fish = fish[fish$agec == fish$Age, ]  
  
# Package the data for stan
  vb_data = list(
    length = fish$Length,
    age = fish$Age,
    obs = seq(1, nrow(fish), 1),
    nobs = nrow(fish),
    group = as.numeric(as.factor(fish$cohort)),
    ngroup = length(unique(fish$cohort)),
    x = as.vector(scale(fish$Year)),
    p_linf = 0,
    p_linf_sd = 1,
    p_k = 0,
    p_k_sd = 1,  
    p_t0 = 0,
    p_t0sd = 1,      
    p_sigma = 10,
    p_tau = 1,
    p_omega = 4,
    p_b_linf_sd = 1,
    p_b_k_sd = 1
  )

# Fit the model with stan
  fit <- stan(file = 'models/vonbert_hmv_cov.stan',
              data = vb_data,
              chains = 3,
              iter = 20000,
              warmup = 15000,
              control = list(
                adapt_delta = 0.90,
                max_treedepth = 10
              )
              )

# Print model summary
  print(fit, digits=3)

# Extract parameters
  pars = extract(fit)  
  
  quantile(pars$ba_linf, c(0.025,0.975))
  quantile(pars$ba_k, c(0.025,0.975))
  
# Leave-one-out cross validation  
  loo(pars$log_lik,
      r_eff=relative_eff(pars$log_lik,
                         chain_id = rep(1:3, 5000))) 

  # loo-ic table without back-calcs
  # ----------------------------------
  # Model    Estimate   SE Linf/K Rank
  # ----------------------------------
  # year       4872.9 39.2    -/+    4 
  # ha         4872.7 40.2    +/-    3 
  # nha        4866.6 40.2    -/ns   2 
  # kg_ha      4867.0 40.1    -/ns   1
  # ----------------------------------
  
# Save result to a file
  save(fit, file='results/vonbert_cohort-x-year.rda')
 
