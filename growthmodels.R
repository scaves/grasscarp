# Front-end needs ---------------------------------------------------------
# Load necessary packages
  library(R2jags)
  library(FSA)

# Clear the global environment
  rm(list=ls())

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
  
  
# Data manipulation -------------------------------------------------------
# Data read
  fish = read.csv('grasscarplengths.csv')  
  
# Drop missing data   
  fish = fish[!is.na(fish$Age) & !is.na(fish$Length), ]
head(fish)
tail(fish)  
nrow(fish)

names(fish)[2] <- 'yearc'


maxes <- aggregate(Age~fishID, fish, max)
names(maxes)[2] <- 'agec'
fish <- merge(fish, maxes)

fish$Year <- fish$yearc-(fish$agec-fish$Age)


# basic von bert ------------------------------------------------------------------


# lets makle a von bert model, and then let it converge for us
vb = nls(
  formula = Length~Linf*(1-exp(-k*(Age-t0))),
  data = fish,
  start = list(Linf = 1250, k = 0.2, t0 = 0), 
  trace = TRUE
)
# now we summarize the model
summary(vb)
# dont report significances, just means and std errors (which are extremely tight
# due to the shear size of the dataset)

# max length of a fish
max(fish$Length)
# make a table of the results, the model coeeficients
res = data.frame(summary(vb)$coefficients)
res
# make a sequence from min to max age, make sure to get rid of the na in dataset
ages = seq(from = min(fish$Age, na.rm = TRUE), 
           to = max(fish$Age, na.rm = TRUE), 
           by = 1)
ages

# make some predictions, using brackets to get values from object res, and plug 
# in new ages just made in the sequence, remake the von bert

meanPred = res[1,1]*(1-exp(-res[2,1]*(ages-res[3,1])))
meanPred

uPred = (res[1,1]+res[1,2]*1.96)*(1-exp(-(res[2,1]+res[2,2]*1.96)*(ages-(res[3,1]+res[3,2]*1.96))))
lPred = (res[1,1]-res[1,2]*1.96)*(1-exp(-(res[2,1]-res[2,2]*1.96)*(ages-(res[3,1]-res[3,2]*1.96))))


# quick check of the values, plot them
plot(x = fish$Age,
     y = fish$Length,
     pch = 21,
     bg = rgb(0.5, 0.5, 0.5, 0.15),
     col = rgb( 0.5, 0.5, 0.5, 0.15),
     xlim = c(0,23),
     ylim = c(0,1400),
     yaxt = "n",
     xlab = 'Age (Years)',
     ylab = 'Length (mm)',
     cex.axis = 1.10,
     cex.lab = 1.15
)
axis(side = 2, las = 2)

lines(x = ages, y = meanPred, lty = 1, col = 'red', lwd = 2) 

# add the confidence intervals
lines(x = ages, y = lPred, lty = 2, col = 'blue', lwd = 2) 
lines(x = ages, y = uPred, lty = 2, col = 'blue', lwd = 2) 
# next week we will bootstrap for confidence intervals, to account for more noise




# bootstrap fitting ------------------------------------------------------------------

# need to make a blank list in order to fill it with the loop
nruns = 1000
out = vector(mode = 'list', length = nruns)

# make the loop
for(i in 1:nruns){
  # make a subset of the data to test each time
  subdat = fish[sample(nrow(fish), 100, replace = F), ]
  
  # remaking the model, just to have it in this section
  
  
  vb = nls(
    formula = Length~Linf*(1-exp(-k*(Age-t0))),
    data = subdat,
    start = list(Linf = 1250, k = 0.2, t0 = 0), 
    trace = FALSE # dont need to see what is happening everytime like above
  )
  
  # make a results table for the output
  out[[i]] = data.frame(summary(vb)$coefficients)[,1]
  
}
# we made the ages sequence earlier, print it to refamiliarize yourself with it

ages

# plot the raw data
plot(x = fish$Age,
     y = fish$Length,
     pch = 21,
     bg = rgb(0.5, 0.5, 0.5, 0.15),
     col = rgb( 0.5, 0.5, 0.5, 0.15),
     xlim = c(0,23),
     ylim = c(0,1400),
     yaxt = "n",
     xlab = 'Age (Years)',
     ylab = 'Length (mm)',
     cex.axis = 1.10,
     cex.lab = 1.15
)
axis(side = 2, las = 2)

# make a second results table, use the function rbind (row bind) to the vector out
res2 = do.call(rbind, out)
res2
# this just takes the vector a puts it all together

# lets make a loop to plot the lines as well so we dont have to type it out 100 times
for(i in 1:nrow(res2)){
  Pred = res2[i,1]*(1-exp(-res2[i,2]*(ages-res2[i,3])))
  
  lines(ages, Pred, col = rgb(0.25, 0.25, 0.25, 0.05))
  
}


# VBGM without site effects, between year captured -----------------------------------
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
  # writeLines(modelString, con='vb_model.txt') 
  # take this line out and dont write a text file so we arent saving to the drive

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
  ni <- 2500  # Number of draws from posterior (for each chain)
  nt <- 500       # Thinning rate
  nb <- 150  # Number of draws to discard as burn-in
  nc <- 3          # Number of chains

# Call jags and run the model
  vb_mod <- jags(data=vb_data, inits=inits, params, textConnection(modelString),
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
    working.directory = getwd())
# make sure we replace the text file with a new function definging the model above

  
# Results of year capture model -----  
# Print a summary of the model
  print(vb_mod)
  
# Get posterior distributions for parameter estimates
  k = vb_mod$BUGSoutput$sims.list$K
  t0 = vb_mod$BUGSoutput$sims.list$to
  linf = vb_mod$BUGSoutput$sims.list$Linf
  


# visualization of the data (year capture) -------------------------------------------------------

boxplot(k, col = 'gray87', xlab = 'year', ylab = 'k') 
  
boxplot(linf, col = 'gray87', xlab = 'year', ylab = 'linf') 
  
# show the different factor levels within the data
fish$YearF = as.numeric(as.factor(fish$Year)) 

# make a sequences of ages for the plots
ages = seq(1, max(fish$Age), 1)
# take the posterior distributions and create the vbgm for each year
first = mean(linf[,1]) * (1-exp(-mean(k[,1])*(ages-mean(t0[,1]))))
second = mean(linf[,2]) * (1-exp(-mean(k[,2])*(ages-mean(t0[,2]))))
third = mean(linf[,3]) * (1-exp(-mean(k[,3])*(ages-mean(t0[,3]))))
fourth = mean(linf[,4]) * (1-exp(-mean(k[,4])*(ages-mean(t0[,4]))))
fifth = mean(linf[,5]) * (1-exp(-mean(k[,5])*(ages-mean(t0[,5]))))
  
  
# now plot the results
# Make the posterior predictive plot, just the raw data, no axes
      par(mar=c(4,4,1,1))
      plot(fish$Age, fish$Length, 
           ylim=c(0, 1500),
           yaxt='n',
           xlab='',
           ylab='',
           xlim=c(0, 25),
           axes = FALSE,
           pch = 21, bg='gray87', col='gray87', main='')
# add the growth curves  
lines(x = ages, y = first, col = 'black')
lines(x = ages, y = second, col = 'red')
lines(x = ages, y = third, col = 'blue')
lines(x = ages, y = fourth, col = 'green')
lines(x = ages, y = fifth, col = 'yellow')

# add in the axes and labels
     axis(1, pos=0)
     axis(2, pos=0, las=2)  
     mtext(expression(paste('Age (years)')),
           side=1, line=2.5)
     mtext('Total length (mm)', side=2, line=2.5) 
legend('bottomright', inset = 0.05,legend=c("2006", "2007", '2009', '2010', '2017'),
       col=c("black", "red", 'blue', 'green', 'yellow'), lty = 1, title = 'Year Capture',
       box.lty = 0)

# figure out the max age of fish in each given year
max(fish$Age[fish$Year == 2006])
max(fish$Age[fish$Year == 2007])
max(fish$Age[fish$Year == 2009])
max(fish$Age[fish$Year == 2010])
max(fish$Age[fish$Year == 2017])

# make the upper confidence interval
firstup = up(linf[,1]) * (1-exp(-up(k[,1])*(ages-up(t0[,1]))))
secondup = up(linf[,2]) * (1-exp(-up(k[,2])*(ages-up(t0[,2]))))
thirdup = up(linf[,3]) * (1-exp(-up(k[,3])*(ages-up(t0[,3]))))
fourthup = up(linf[,4]) * (1-exp(-up(k[,4])*(ages-up(t0[,4]))))
fifthup = up(linf[,5]) * (1-exp(-up(k[,5])*(ages-up(t0[,5]))))
# make the lower confidence interval
firstlow = low(linf[,1]) * (1-exp(-low(k[,1])*(ages-low(t0[,1]))))
secondlow = low(linf[,2]) * (1-exp(-low(k[,2])*(ages-low(t0[,2]))))
thirdlow = low(linf[,3]) * (1-exp(-low(k[,3])*(ages-low(t0[,3]))))
fourthlow = low(linf[,4]) * (1-exp(-low(k[,4])*(ages-low(t0[,4]))))
fifthlow = low(linf[,5]) * (1-exp(-low(k[,5])*(ages-low(t0[,5]))))


# Biomass data read and manipulation -----------------------------------------------------------------

# start by merging all the data together
hydrilla = read.csv('gcStockingAndHydrilla.csv')

hydrilla
names(hydrilla)
names(fish)
colnames(hydrilla)[1] = 'Year'
names(hydrilla)

biomass = merge(x = hydrilla, y = fish, by = 'Year')
head(biomass)
tail(biomass)


# now we have to add biomass as a continuous variable of hydrilla to the model


# Biomass model -----------------------------------------------------------

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
ni1 <- 30000  # Number of draws from posterior (for each chain)
nt1 <- 200       # Thinning rate
nb1 <- 10000  # Number of draws to discard as burn-in
nc1 <- 3          # Number of chains

# Call jags and run the model
vb_mod_cont <- jags(data=vb_data_cont, inits=inits1, params1, textConnection(modelString1),
               n.chains = nc1, n.thin = nt1, n.iter = ni1, n.burnin = nb1,
               working.directory = getwd())
# make sure we replace the text file with a new function definging the model above
vb_mod_cont

save(vb_mod_cont, file='covRes.rda')

# data visualization biomass, still a work in progress -----  


# Print a summary of the model
print(vb_mod_cont)

beta0_l = vb_mod_cont$BUGSoutput$sims.list$beta0_l
betah_l = vb_mod_cont$BUGSoutput$sims.list$betah_l
beta0_k = vb_mod_cont$BUGSoutput$sims.list$beta0_k
betah_k = vb_mod_cont$BUGSoutput$sims.list$betah_k

# exp(lLinf)


# vb_mod_cont$BUGSoutput$sims.list



# plotting code that I haven't looked at yet becasue the model isn't working
# next step is to tear this apart and determine best way to plot continuous 
# variable (hydrilla biomass) model



# since all data is on range from -2 to 2
newBiomass = seq(-2, 2, 0.1)    

# make an empty matrix to then fill with data later
preds = matrix(NA, nrow = length(beta0), ncol = length(newBiomass))

for(i in 1:nrow(preds)){
  for(t in 1:length(newBiomass)){
    preds[i, t] = (beta0[i] + betah[i]*newBiomass[t])
  }
}
preds = exp(preds)



par(mar = c(5,5,1,1))
plot(x = newBiomass, y = preds[1, ], type = 'l',
    # col = rgb(0.4, 0.4, 0.4, 0.05),
     xlab = "Biomass", ylab = "k",
     yaxt = 'n', ylim = c(0,1))  
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



