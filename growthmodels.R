# Front-end needs ----
# Load necessary packages
  library(R2jags)
  library(FSA)

# Create function to invert logit link function
  inv.logit = function(x){
    exp(x)/(1+exp(x))
  }

# Make a function to get lower 95% credible limit with short name
  low = function(x){
    quantile(x, probs=c(0.025))
  }

# Make a function to get upper 95% credible limit with short name
  up = function(x){
    quantile(x, probs=c(0.975))
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
biomass = merge(x = hydrilla, y = fish, by = 'Year')

# VBGM by year captured -----
# . Model string -----
  modelString = "
    model{
      # Likelihood
        for(i in 1:N){
          Y[i] ~ dnorm(L[i], tau[Ti[i]])
          L[i] <- Linf[group[i]]*(1-exp(-K[group[i]]*(Ti[i]-to[group[i]])))
        }

      # Priors on VBGM parameters, estimating the parameters for each group (each year)
        for(j in 1:ngroups){ 
        Linf[j] ~ dlnorm(10, 0.001)
        logit(K[j]) <- k[j]
        k[j] ~ dnorm(0, 0.001)
        to[j] ~ dunif(-10, 0)
        }

      # Precision for length at age distribution
        for(t in 1:Tmax){
          tau[t] ~ dgamma(0.001, 0.0001)
        }
    }"

# . JAGS settings -----
# Package the data for JAGS, adding new data for groups and ngroups
  vb_data = list(
    Y = fish$Length,
    Ti = fish$Age,
    N = nrow(fish),
    Tmax = max(fish$Age),
    group = as.numeric(as.factor(fish$yearc)),
    ngroups = length(unique(fish$yearc))
  )

# Parameters monitored
  params = c('Linf', 'K', 'to')

# Initial values for parameters
  inits <- function(){
    list(
      Linf = runif(length(unique(fish$yearc)), 1, 10),
      k = rnorm(length(unique(fish$yearc)), 0, 1),
      to = runif(length(unique(fish$yearc)), -10, 0),
      tau = rgamma(max(fish$Age), .01, 1)
    )
  }

# MCMC settings
  # need to change iterations, thinning rate and burnins to get the model to 
  # converge, it is close and I am getting good values, dave used 500000 iterations,
  # and thinning rate of 100
  ni <- 250000  # Number of draws from posterior (for each chain)
  nt <- 200     # Thinning rate
  nb <- 100000  # Number of draws to discard as burn-in
  nc <- 3       # Number of chains

# . Model calibration -----
# Call jags and run the model
  vb_mod <- jags(data=vb_data, inits=inits, params, textConnection(modelString),
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
    working.directory = getwd())
  
# Save the model object to a file
  save(vb_mod, file='yearRes.rda')
  
  
# . Results of year capture model ----- 
# Load data
  load("yearRes.rda")  
  
# Print a summary of the model
  print(vb_mod)
  
# Get posterior distributions for parameter estimates
  k = vb_mod$BUGSoutput$sims.list$K
  t0 = vb_mod$BUGSoutput$sims.list$to
  linf = vb_mod$BUGSoutput$sims.list$Linf
  
# .. Data visualization by year -----
# ... Boxplot of k by year of capture -----
# Set graphical parameters
par(mar=c(5,5,1,1))
  
# Make initial boxplot to get 
# boxplot object, saved to h
h=boxplot(k, plot=FALSE)

# Replace whisker length with the
# 95% Confidence intervals
h$stats[c(1,5), ] = apply(k, 2, quantile, probs=c(0.025, 0.975))[c(1,2), ]

# Re-plot the boxplot with the new whiskers 
# and add graphical parameters to make the
# plot prettier
bxp(h, boxfill='gray87', outline=FALSE, ylim=c(0,.25), xaxt='n',
    yaxt='n', plot=FALSE, boxwex=.5,
    pars = list(staplewex=0, whisklty=1, whiskcol="gray40", whisklwd=2,
                boxcol='gray40', boxfill='gray87', medcol='gray40')
    )
box()
# Add axes and axis labels
axis(1, at=seq(1,5,1), sort(unique(fish$yearc)))
axis(2, las=2, cex.axis=1.10)
mtext(side = 1, 'Year of Collection', line=3.5, cex=1.15)
mtext(side = 2, 'K', line=3.5, cex=1.15)

# ... Boxplot of linf by year of capture -----
# Set graphical parameters
par(mar=c(5,5,1,1))
# Make initial boxplot to get 
# boxplot object, saved to h
h=boxplot(linf, plot=FALSE)
# Replace whisker length with the
# 95% Confidence intervals
h$stats[c(1,5), ] = apply(linf, 2, quantile, probs=c(0.025, 0.975))[c(1,2), ]
# Re-plot the boxplot with the new whiskers 
# and add graphical parameters to make the
# plot prettier
bxp(h, boxfill='gray87', outline=FALSE, ylim=c(0,3000), xaxt='n',
          yaxt='n', plot=FALSE, boxwex=.5,
    pars = list(staplewex=0, whisklty=1, whiskcol="gray40", whisklwd=2,
                boxcol='gray40', boxfill='gray87', medcol='gray40')
    )
box()
axis(1, at=seq(1,5,1), sort(unique(fish$yearc)))
axis(2, las=2, cex.axis=1.10)
mtext(side = 1, 'Year of Collection', line=3.5, cex=1.15)
mtext(side = 2, expression(paste('L'[infinity])), line=3.5, cex=1.15)
  

# ... Growth curves by year -----
# Make a sequences of ages for the plots
ages = seq(1, max(fish$Age), 1)

# Posterior predictive vbgm for each year
first = mean(linf[,1]) * (1-exp(-mean(k[,1])*(ages-mean(t0[,1]))))
second = mean(linf[,2]) * (1-exp(-mean(k[,2])*(ages-mean(t0[,2]))))
third = mean(linf[,3]) * (1-exp(-mean(k[,3])*(ages-mean(t0[,3]))))
fourth = mean(linf[,4]) * (1-exp(-mean(k[,4])*(ages-mean(t0[,4]))))
fifth = mean(linf[,5]) * (1-exp(-mean(k[,5])*(ages-mean(t0[,5]))))
  
# Plot the results
# First, the raw data with no axes
par(mar=c(4,4,1,1))
plot(fish$Age, fish$Length, 
     ylim=c(0, 1500),
     yaxt='n',
     xlab='',
     ylab='',
     xlim=c(0, 25),
     axes = FALSE,
     pch = 21,
     bg=c("black", "red", 'blue', 'green', 'orange')[as.factor(fish$yearc)],
     col=c("black", "red", 'blue', 'green', 'orange')[as.factor(fish$yearc)],
     main='')

# Add the growth curves  
lines(x = ages, y = first, col = 'black', lwd = 2)
lines(x = ages, y = second, col = 'red', lwd = 2)
lines(x = ages, y = third, col = 'blue', lwd = 2)
lines(x = ages, y = fourth, col = 'green', lwd = 2)
lines(x = ages, y = fifth, col = 'orange', lwd = 2)

# Add axes and labels
axis(1, pos=0)
axis(2, pos=0, las=2)  
mtext(expression(paste('Age (years)')), side=1, line=2.5)
mtext('Total length (mm)', side=2, line=2.5) 
 
# Add legend
legend('bottomright', inset = 0.05,
       legend=c("2006", "2007", '2009', '2010', '2017'),
       col=c("black", "red", 'blue', 'green', 'orange'),
       lty = 1, title = 'Year of Capture', box.lty = 0, lwd = 2)


# Hydrilla model ----
# . Model string ----
modelString1 = "
    model{
      # Likelihood
          for(i in 1:N){
            Y[i] ~ dnorm(L[i], tau[Ti[i]])
            L[i] <- lLinf[i]*(1-exp(-K[i]*(Ti[i]-to)))
            
            log(K[i]) <- beta0_k + betah_k*ha[i]
            log(lLinf[i]) <- beta0_l + betah_l*ha[i]
          
          }

      # Priors on VBGM parameters, estimating the parameters for each group 
        betah_k ~ dnorm(0, 0.0001)
        beta0_k ~ dnorm(0, 0.0001)
        beta0_l ~ dnorm(6.5, 0.0001)
        betah_l ~ dnorm(0, 0.0001)
        to ~ dnorm(0, 0.0001)

      # Precision for length at age distribution
        for(t in 1:Tmax){
          tau[t] ~ dgamma(0.001, 0.0001)
        }
     }"


# . JAGS settings ----
# Package the data for JAGS, adding new data for groups and ngroups
vb_data_cont = list(
  Y = biomass$Length,
  ha = as.vector(scale(biomass$ha)),
  Ti = biomass$Age,
  N = nrow(biomass),
  Tmax = max(biomass$Age)
)

# Parameters monitored
params1 = c('beta0_k', 'betah_k', 'betah_l', 'beta0_l', 'to')

# Initial values for parameters
inits1 <- function(){
  list(
    beta0_k = rnorm(1, 0, 1),
    betah_k = rnorm(1, 0, 1),
    betah_l = rnorm(1, 0, 1),
    beta0_l = rnorm(1, 6.5, .01),
    to = runif(1, -10, 0),
    tau = rgamma(max(fish$Age), .01, 1)
  )
}

# MCMC settings
# need to change iterations, thinning rate and burnins to get the model to 
# converge, it is close and I am getting good values, dave used 500000 iterations,
# and thinning rate of 100
ni1 <- 250000 # Number of draws from posterior (for each chain)
nt1 <- 200    # Thinning rate
nb1 <- 10000  # Number of draws to discard as burn-in
nc1 <- 3      # Number of chains

# .. Model calibration ----
# Call jags and run the model
vb_mod_cont <- jags(data=vb_data_cont, inits=inits1, params1, textConnection(modelString1),
               n.chains = nc1, n.thin = nt1, n.iter = ni1, n.burnin = nb1,
               working.directory = getwd())

# Save the model results
save(vb_mod_cont, file='covRes.rda')
# . Results ----
# .. Data visualization for hydrilla model ----  
# Load data
load("covRes.rda")
  
# Print a summary of the model
print(vb_mod_cont)

beta0_l = vb_mod_cont$BUGSoutput$sims.list$beta0_l
betah_l = vb_mod_cont$BUGSoutput$sims.list$betah_l
beta0_k = vb_mod_cont$BUGSoutput$sims.list$beta0_k
betah_k = vb_mod_cont$BUGSoutput$sims.list$betah_k

# .. Plotting code for effect of hydrilla -----
# ... Hdrilla vs L-infinity -----
# since all data is on range from -2 to 2
newBiomass = seq(-2, 2, 0.1)    

# make an empty matrix to then fill with data later
preds = matrix(NA, nrow = length(beta0_l), ncol = length(newBiomass))
# Fill the matrix with predicted L-inf
# based on the linear predictor that
# uses hydrilla biomass
for(i in 1:nrow(preds)){
  for(t in 1:length(newBiomass)){
    preds[i, t] = (beta0_l[i] + betah_l[i]*newBiomass[t])
  }
}
preds = exp(preds)


# Plot the first prediction to get a plotting window set up
par(mar = c(5,5,1,1))
plot(x = newBiomass, y = preds[1, ], type = 'l',
     col = rgb(0.4, 0.4, 0.4, 0.05),
     xlab = "Standardized Surface Hectares", ylab = expression(paste('L'[infinity])),
     yaxt = 'n', ylim = c(1000,1200))  
axis(2, las = 2)
# now we add the loop for the rest of the data
for (i in 2:nrow(preds)){
  lines(x = newBiomass, y = preds[i, ],
        col = rgb(0.4, 0.4, 0.4, 0.05))
}

# get the mean and 95% CRIs, defined the upper and lower CRI functions above
lines(x = newBiomass, y = apply(preds, MARGIN = 2, FUN = mean),
      col = 'blue', lwd = 2)
lines(x = newBiomass, y = apply(preds, MARGIN = 2, FUN = low),
      col = 'red', lwd = 2, lty = 2)
lines(x = newBiomass, y = apply(preds, MARGIN = 2, FUN = up),
      col = 'red', lwd = 2, lty = 2)
box()

# ... Hydrilla vs k -----
# since all data is on range from -2 to 2
newBiomass = seq(-2, 2, 0.1)    

# make an empty matrix to then fill with data later
preds = matrix(NA, nrow = length(beta0_k), ncol = length(newBiomass))
# Fill the matrix with predicted L-inf
# based on the linear predictor that
# uses hydrilla biomass
for(i in 1:nrow(preds)){
  for(t in 1:length(newBiomass)){
    preds[i, t] = (beta0_k[i] + betah_k[i]*newBiomass[t])
  }
}
preds = exp(preds)

# Plot the first prediction to get a plotting window set up
par(mar = c(5,5,1,1))
plot(x = newBiomass, y = preds[1, ], type = 'l',
     col = rgb(0.4, 0.4, 0.4, 0.05),
     xlab = "Standardized Surface Hectares", ylab = 'K',
     yaxt = 'n', ylim = c(0.1, 0.3))  
axis(2, las = 2)
# now we add the loop for the rest of the data
for (i in 2:nrow(preds)){
  lines(x = newBiomass, y = preds[i, ],
        col = rgb(0.4, 0.4, 0.4, 0.05))
}

# get the mean and 95% CRIs, defined the upper and lower CRI functions above
lines(x = newBiomass, y = apply(preds, MARGIN = 2, FUN = mean),
      col = 'blue', lwd = 2)
lines(x = newBiomass, y = apply(preds, MARGIN = 2, FUN = low),
      col = 'red', lwd = 2, lty = 2)
lines(x = newBiomass, y = apply(preds, MARGIN = 2, FUN = up),
      col = 'red', lwd = 2, lty = 2)
box()



# ... Hydrilla vs M -----
# since all data is on range from -2 to 2
newBiomass = seq(-2, 2, 0.1)    

# make an empty matrix to then fill with data later
preds = matrix(NA, nrow = length(beta0_k), ncol = length(newBiomass))
# Fill the matrix with predicted L-inf
# based on the linear predictor that
# uses hydrilla biomass
for(i in 1:nrow(preds)){
  for(t in 1:length(newBiomass)){
    preds[i, t] = (beta0_k[i] + betah_k[i]*newBiomass[t])
  }
}
preds = exp(preds)
preds <- apply(preds, 2, FUN='*', 1.5)

# Plot the first prediction to get a plotting window set up
par(mar = c(5,5,1,1))
plot(x = newBiomass, y = preds[1, ], type = 'l',
     col = rgb(0.4, 0.4, 0.4, 0.05),
     xlab = "Standardized Surface Hectares", ylab = 'M',
     yaxt = 'n', ylim = c(0.1, 0.3))  
axis(2, las = 2)
# now we add the loop for the rest of the data
for (i in 2:nrow(preds)){
  lines(x = newBiomass, y = preds[i, ],
        col = rgb(0.4, 0.4, 0.4, 0.05))
}

# get the mean and 95% CRIs, defined the upper and lower CRI functions above
lines(x = newBiomass, y = apply(preds, MARGIN = 2, FUN = mean),
      col = 'blue', lwd = 2)
lines(x = newBiomass, y = apply(preds, MARGIN = 2, FUN = low),
      col = 'red', lwd = 2, lty = 2)
lines(x = newBiomass, y = apply(preds, MARGIN = 2, FUN = up),
      col = 'red', lwd = 2, lty = 2)
box()