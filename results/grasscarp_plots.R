# Front-end ----
# Package loads
  library(rstan)


# . Function definition ----
# Make a function to get lower 95% credible limit with short name
  low = function(x){
    quantile(x, probs=c(0.025))
  }

# Make a function to get upper 95% credible limit with short name
  up = function(x){
    quantile(x, probs=c(0.975))
  }
  
# Scale new values of variable
# based on a sample
  nscale <- function(xn, xo){
    (xn - mean(xo) ) / sd(xo)
  }

# Unscale
  unscale <- function(x, y){
    x * sd(y) + mean(y)
  }    
  
# Make function to set negative values
# to zero
  zeroes <- function(x){
    x[x<0] <- 0
    return(x)
  }  
  
# Data manipulation ----
# . Fish data ----
# Read in length-at-age data
  fish = read.csv('data/grasscarplengths.csv')  
  
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
  fish = merge(x = hydrilla, y = fish, by = 'Year')

  
# Data visualization ----
# . Year model ----
# .. Load results ----
# Load the model results
  load('results/yearmod-result.rda')

# Get parameters from the model
  parms <- extract(fit)
   
# Get posterior distributions for parameter estimates
  k = parms$K
  t0 = parms$t0
  linf = parms$Linf  
  
# .. Boxplots of posteriors by year -----
# Set up an image file
tiff(filename='results/Figure1.tiff',
     width = 1500,
     height = 2000,
     res = 500,
     pointsize = 8,
     units = 'px'
     )

# Set graphical parameters  
par(mfrow=c(3,1), oma=c(4,5,1,1), mar=c(1,1,1,1))

# Linf by year of capture 
  # Make initial boxplot to get 
  # boxplot object, saved to h
    h <- boxplot(linf, plot=FALSE)
  # Replace whisker length with the
  # 95% Confidence intervals
    h$stats[c(1,5), ] <- 
      apply(linf, 2, quantile, probs=c(0.025, 0.975))[c(1,2), ]
  # Re-plot the boxplot with the new whiskers 
  # and add graphical parameters to make the
  # plot prettier
    bxp(h, boxfill='gray87',
        outline=FALSE,
        ylim=c(0,3000),
        xaxt='n', yaxt='n',
        plot=FALSE,
        boxwex=.5, 
        pars = list(
          staplewex=0,
          staplecol='gray40',
          whisklty=1,
          whiskcol="gray40",
          whisklwd=1,
          boxcol='gray40',
          boxfill='gray87',
          medcol='gray40',
          medlwd=1
          )
        )
    box()
  axis(1, labels = FALSE, tick = TRUE)
  axis(2, las = 2)
  mtext(side = 2, expression(paste('L'[infinity])), line=3.5)

  # K by year of capture
  # Make initial boxplot to get 
  # boxplot object, saved to h
    h=boxplot(k, plot=FALSE)
  # Replace whisker length with
  # 95% Confidence intervals
    h$stats[c(1,5), ] <- 
      apply(k, 2, quantile, probs=c(0.025, 0.975))[c(1,2), ]
  # Re-plot the boxplot with the new whiskers 
  # and add graphical parameters to make the
  # plot prettier
    bxp(h, boxfill='gray87',
        outline=FALSE,
        ylim=c(0,.25),
        xaxt='n', yaxt='n',
        plot=FALSE, 
        pars = list(
          staplewex=0,
          staplecol='gray40',
          whisklty=1, 
          whiskcol="gray40", 
          whisklwd=1,
          boxcol='gray40', 
          boxfill='gray87', 
          boxwex=.5, 
          medcol='gray40',
          medlwd=1
          )
        )
    box()
  # Add axes and axis labels
    axis(1, labels=FALSE, tick=TRUE)
    axis(2, las=2, cex.axis=1.10)
    mtext(side = 2, 'K', line=3.5)
    
# t0 by year of capture 
  # Make initial boxplot to get 
  # boxplot object, saved to h
    h=boxplot(t0, plot=FALSE)
  # Replace whisker length with the
  # 95% Confidence intervals
    h$stats[c(1,5), ] <- 
      apply(t0, 2, quantile, probs=c(0.025, 0.975))[c(1,2), ]
  # Re-plot the boxplot with the new whiskers 
  # and add graphical parameters to make the
  # plot prettier
    bxp(h, boxfill='gray87',
        outline=FALSE,
        ylim=c(-3,0),
        xaxt='n', yaxt='n',
        plot=FALSE,
        boxwex=.5, 
        pars = list(
          staplewex=0,
          staplecol='gray40',
          whisklty=1,
          whiskcol="gray40",
          whisklwd=1,
          boxcol='gray40',
          boxfill='gray87',
          medcol='gray40',
          medlwd=1
          )
        )
    box()
  axis(1, at=seq(1,5,1), sort(unique(fish$yearc)))
  axis(2, las=2)
  mtext(side = 1, 'Year of Collection', line=3.5)
  mtext(side = 2, expression(paste('t'['0'])), line=3.5)

dev.off()    
  
# .. Growth curves by year -----
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
     bg='gray87',
     col='gray87',
     main='')

# Add the growth curves  
lines(x = ages, y = first, col = 'black', lwd = 2, lty = 1)
lines(x = ages, y = second, col = 'black', lwd = 2, lty = 2)
lines(x = ages, y = third, col = 'black', lwd = 2, lty = 3)
lines(x = ages, y = fourth, col = 'black', lwd = 2, lty = 4)
lines(x = ages, y = fifth, col = 'black', lwd = 2, lty = 5)

# Add axes and labels
axis(1, pos=0)
axis(2, pos=0, las=2)  
mtext(expression(paste('Age (years)')), side=1, line=2.5)
mtext('Total length (mm)', side=2, line=2.5) 
 
# Add legend
legend('bottomright', inset = 0.05,
       legend=c("2006", "2007", '2009', '2010', '2017'),
       col='black',
       lty = c(1, 2, 3, 4, 5), title = 'Year of Capture', box.lty = 0, lwd = 2)


# . Hydrilla MV model -----
# .. Data ----
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

# .. Load the results ----
# Load the model fit object called 'fit'
load("results/vonbert_hydrilla_mv.rda") 

# Extract the model parameters
pars = extract(fit)
  
# .. Set up an image file ----
tiff(filename='results/Figure2.tiff',
     width = 1500,
     height = 2000,
     res = 500,
     pointsize = 8,
     units = 'px'
     )

# Set graphical parameters  
par(mfrow=c(3,1), oma=c(4,5,1,1), mar=c(1,1,1,1))

# .. Hdrilla vs L-infinity -----
# Make a sequence of new biomasses
newBiomass = seq(min(vb_data$hydrilla), max(vb_data$hydrilla), 0.1)

# Make an empty matrix to then fill with data later
preds = matrix(NA, nrow = length(pars$b0_linf), ncol = length(newBiomass))

# Fill the matrix with predicted L-inf
# based on the linear predictor that
# uses hydrilla biomass
for(i in 1:nrow(preds)){
  for(t in 1:length(newBiomass)){
    preds[i, t] = exp(pars$mu_beta_cor[i,1] + pars$b0_linf[i] + pars$bh_linf[i]*newBiomass[t])
  }
}

# Plot the first prediction to get a plotting window set up
plot(x = newBiomass, y = preds[1, ], type = 'l',
     col = rgb(0.6, 0.6, 0.6, 0.05),
     xlab = "", 
     xaxt = 'n',
     xlim = c(min(vb_data$hydrilla), max(vb_data$hydrilla)),
     ylab = '',
     yaxt = 'n',
     ylim = c(1000,2000)
     )  
mtext(side=2, expression(paste('L'[infinity])), line=3.5)
axis(2, las = 2)
axis(1, at=nscale(seq(0,1400,200), fish$ha), labels=FALSE)

# Calculate the mean and 95% CRIs for posterior predictions
muPred = apply(preds, 2, mean)
lowPred = apply(preds, 2, low)
lowPred = zeroes(lowPred)
upPred = apply(preds, 2, up)

polygon(c(newBiomass, rev(newBiomass)), c(lowPred, rev(upPred)),
        col=rgb(.8, .8, .8, 0.5),
        border=NA
        )

# get the mean and 95% CRIs, defined the upper and lower CRI functions above
lines(x = newBiomass, y = muPred, col = 'black', lwd = .5)
#lines(x = newBiomass, y = lowPred, col = 'black', lwd = 1, lty = 2)
#lines(x = newBiomass, y = upPred, col = 'black', lwd = 1, lty = 2)


# .. Hydrilla vs k -----
# Make a sequence of new biomasses
newBiomass = seq(min(vb_data$hydrilla), max(vb_data$hydrilla), 0.1)

# Make an empty matrix to then fill with data later
preds = matrix(NA, nrow = length(pars$b0_k), ncol = length(newBiomass))
# Fill the matrix with predicted L-inf
# based on the linear predictor that
# uses hydrilla biomass
for(i in 1:nrow(preds)){
  for(t in 1:length(newBiomass)){
    preds[i, t] = exp(pars$mu_beta_cor[i,2] + pars$b0_k[i] + pars$bh_k[i]*newBiomass[t])
  }
}

# Plot the first prediction to get a plotting window set up
plot(x = newBiomass, y = preds[1, ], type = 'l',
     col = rgb(0.6, 0.6, 0.6, 0.05),
     xaxt='n',
     xlim = c(min(vb_data$hydrilla), max(vb_data$hydrilla)),
     xlab = "",
     ylab = '',
     yaxt = 'n',
     ylim = c(0, 0.5)
     )
axis(1, at=nscale(seq(0,1400,200), fish$ha), labels=FALSE)
mtext(side=2, 'K', line=3.5)
axis(2, las = 2)

# Calculate the mean and 95% CRIs for posterior predictions
muPred = apply(preds, 2, mean)
lowPred = apply(preds, 2, low)
lowPred = zeroes(lowPred)
upPred = apply(preds, 2, up)

polygon(c(newBiomass, rev(newBiomass)), c(lowPred, rev(upPred)),
        col=rgb(.8, .8, .8, 0.5),
        border=NA
        )

# get the mean and 95% CRIs, defined the upper and lower CRI functions above
lines(x = newBiomass, y = muPred, col = 'black', lwd = .5)
#lines(x = newBiomass, y = lowPred, col = 'black', lwd = 1, lty = 2)
#lines(x = newBiomass, y = upPred, col = 'black', lwd = 1, lty = 2)


# .. Hydrilla vs M -----
# Make an empty matrix to then fill with data later
preds = matrix(NA, nrow = length(pars$b0_linf), ncol = length(newBiomass))
# Fill the matrix with predicted L-inf
# based on the linear predictor that
# uses hydrilla biomass
for(i in 1:nrow(preds)){
  for(t in 1:length(newBiomass)){
    preds[i, t] = exp(pars$mu_beta_cor[i,2] + pars$b0_k[i] + pars$bh_k[i]*newBiomass[t])
  }
}
preds <- apply(preds, 2, FUN='*', 1.5)

# Plot the first prediction to get a plotting window set up
plot(x = newBiomass, y = preds[1, ], type = 'l',
     col = rgb(0.6, 0.6, 0.6, 0.05),
     xlim = c(min(vb_data$hydrilla), max(vb_data$hydrilla)),
     xlab = '', xaxt='n',ylab = '',
     yaxt = 'n', ylim = c(0, .5)) 
mtext(side=1, "Surface hectares of hydrilla", line=3.5)
mtext(side=2, 'M', line=3.5)
axis(1, at=nscale(seq(0,1400,200), fish$ha),
     labels=seq(0,1400,200))
axis(2, at=seq(0,1,.1), labels=seq(0,1,.1), las = 2)

# Calculate the mean and 95% CRIs for posterior predictions
muPred = apply(preds, 2, mean)
lowPred = apply(preds, 2, low)
lowPred = zeroes(lowPred)
upPred = apply(preds, 2, up)

polygon(c(newBiomass, rev(newBiomass)), c(lowPred, rev(upPred)),
        col=rgb(.8, .8, .8, 0.5),
        border=NA
        )

# get the mean and 95% CRIs, defined the upper and lower CRI functions above
lines(x = newBiomass, y = muPred, col = 'black', lwd = .5)
#lines(x = newBiomass, y = lowPred, col = 'black', lwd = 1, lty = 2)
#lines(x = newBiomass, y = upPred, col = 'black', lwd = 1, lty = 2)

dev.off()

# .. Growth curve ----
# Load the results
fish = fish[fish$agec == fish$Age, ]
load("results/vonbert_hydrilla_mv.rda")

# Extract parameter estimates from model fit
pars = extract(fit)

# Offsets
Linf = exp(pars$mu_beta_cor[,1] + pars$b0_linf)
K = exp(pars$mu_beta_cor[,2] + pars$b0_k)
t0 = pars$mu_beta_cor[,3] + pars$b0_t0

# Make a sequence of new ages for which we will predict lengths
Age = seq(0, 23, .1)

# Predict mean length at age for each sample
preds = matrix(0, length(Linf), length(Age))
for(i in 1:length(Linf)){
  for(t in 1:length(Age)){
    preds[i,t] = Linf[i]*(1-exp(-K[i]*(Age[t]-(t0[i]))))
  }
}

# Open file connection
tiff(
  filename = paste0('results/Figure1.tiff'),
  height = 878,
  width = 1314,
  pointsize = 32
  )

# Make the plot
par(mar=c(5,5,1,1))
plot(fish$Age, fish$Length, ylim=c(0, 1400),
     yaxt='n', xlab='Age (years)',
     ylab='Fork length (mm)',
     xlim=c(0,23), axes = FALSE,
     pch = 21, bg='white', col='white',
     main='')

# Calculate the mean and 95% CRIs for posterior predictions
muPred = apply(preds, 2, mean)
lowPred = apply(preds, 2, low)
lowPred = zeroes(lowPred)
upPred = apply(preds, 2, up)

# Add a polygon for 95% CI
polygon(c(Age, rev(Age)), c(lowPred, rev(upPred)),
        col=rgb(.8, .8, .8, 0.5),
        border=rgb(.8, .8, .8, 0.5)
        )
# Add the raw data
points(fish$Age, fish$Length,
       pch = 21,
       col=rgb(.4, .4, .4, 0.9),
       bg=rgb(.4, .4, .4, 0.5)
)

# Plot the mean and 95% CRI for predicted length at each age
lines(Age, muPred, col='blue', lwd=1, lty=1)
lines(Age, upPred, col='red', lwd=1, lty=2)
lines(Age, lowPred, col='red', lwd=1, lty=2)

# Plot axes
axis(1, pos=0, at=seq(0,23,1), labels=seq(0,23,1))
axis(2, pos=0, las=2)

# Close file connection
dev.off()  
