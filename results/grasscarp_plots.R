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
  fish1 = fish
  
# Drop back-calculated lengths
  fish = fish[fish$Age == fish$agec,]
  
  
# Data visualization ----
# . Figure 1. Annual VBGF parameters ----
# .. Load results ----
# Load the model results
  load('results/vonbert_annual.rda')

# Get parameters from the model
  parms <- extract(fit)
   
# Get posterior distributions for parameter estimates
  k = parms$K
  t0 = parms$t0
  linf = parms$Linf  
  
# .. Boxplots of posteriors by year -----
# Set up an image file
tiff(filename='results/Figure1.tiff',
     width = 1300,
     height = 2000,
     res = 400,
     pointsize = 8
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
        ylim=c(1000,1800),
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
  mtext(side = 2, expression(paste(hat(italic('L'))[infinity])), line=3.5, cex=.66)

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
        ylim=c(0,.5),
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
    axis(2, las=2)
    mtext(side = 2, expression(hat(italic(K))), line=3.5, cex=.66)
    
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
        ylim=c(-6,0),
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
  mtext(side = 1, 'Year of Collection', line=3.5, cex=.66)
  mtext(side = 2, expression(paste(hat(italic('t'))['0'])), line=3.5, cex=.66)

dev.off()    
  

# . Figure 2. Growth curves -----
# . Set up an image file ----
tiff(filename='results/Figure2.tiff',
     width = 1300,
     height = 2000,
     res = 400,
     pointsize = 5
     )
# Set graphical parameters
  par(mfrow=c(2, 1), oma=c(4,4,1,1), mar=c(1,1,1,1))
  
# .. Year-of-capture model -----
# Load the model results
  load('results/vonbert_annual.rda')  
  
# Get parameters from the model
  parms <- extract(fit)
   
# Get posterior distributions for parameter estimates
  k = parms$K
  t0 = parms$t0
  linf = parms$Linf  

# Make a sequence of new ages for which we will predict lengths
  ages = seq(0, 23, .1)

# Posterior predictive vbgm for each year
  first = mean(linf[,1]) * (1-exp(-mean(k[,1])*(ages-mean(t0[,1]))))
  second = mean(linf[,2]) * (1-exp(-mean(k[,2])*(ages-mean(t0[,2]))))
  third = mean(linf[,3]) * (1-exp(-mean(k[,3])*(ages-mean(t0[,3]))))
  fourth = mean(linf[,4]) * (1-exp(-mean(k[,4])*(ages-mean(t0[,4]))))
  fifth = mean(linf[,5]) * (1-exp(-mean(k[,5])*(ages-mean(t0[,5]))))
    
# Plot the results
  plot(fish$Age, fish$Length, ylim=c(0, 1400),
       yaxt='n', ylab='', xlab='',
       xlim=c(0,25), axes = FALSE,
       pch = 21,
       col=rgb(.4, .4, .4, 0.5),
       bg=rgb(.4, .4, .4, 0.5),
       main='')

# Add the growth curves  
  lines(x = ages, y = first, col = 'black', lwd = .9, lty = 1)
  lines(x = ages, y = second, col = 'black', lwd = .9, lty = 2)
  lines(x = ages, y = third, col = 'black', lwd = .9, lty = 3)
  lines(x = ages, y = fourth, col = 'black', lwd = .9, lty = 4)
  lines(x = ages, y = fifth, col = 'black', lwd = .9, lty = 5)

# Add axes and labels
  axis(1, pos=0, labels=FALSE)
  axis(2, pos=0, las=2)  
   
# Add legend
  legend('bottomright',
         inset = 0.05,
         legend=c("2006", "2007", '2009', '2010', '2017'),
         col='black',
         lty = c(1, 2, 3, 4, 5),
         title = 'Year of Capture',
         box.lty = 0,
         lwd = .9
         )


# .. Covariate model ----
# Load the results
  load("results/vonbert_hydrilla_mv.rda")

# Extract parameter estimates from model fit
  pars = extract(fit)

# Create random draw from standard normal to 
# include uncertainty related to standardized
# covariate
  covs <- rnorm(length(pars$b0_k), 0, 1)
  
# Offsets and covariate uncertainty
  Linf = exp(pars$Gamma[,1] + pars$b0_linf + pars$bh_linf*covs)
  K = exp(pars$Gamma[,2] + pars$b0_k + pars$bh_k*covs)
  t0 = exp(pars$Gamma[,3] + pars$b0_t0)-10

# Make a sequence of new ages for which we will predict lengths
  Age = seq(0, 23, .1)

# Predict mean length at age for each sample
  preds = matrix(0, length(Linf), length(Age))
  for(i in 1:length(Linf)){
    for(t in 1:length(Age)){
      preds[i,t] = Linf[i]*(1-exp(-K[i]*(Age[t]-(t0[i]))))
    }
  }

# Make the plot
  plot(fish$Age, fish$Length, 
       ylab = '',
       ylim=c(0, 1400),
       yaxt='n', 
       xlab='',
       xlim=c(0,25),
       axes = FALSE,
       pch = 21,
       col=rgb(.4, .4, .4, 0.5),
       bg=rgb(.4, .4, .4, 0.5),
       main='')

# Calculate the mean and 95% CRIs for posterior predictions
  muPred = apply(preds, 2, mean)
  lowPred = apply(preds, 2, low)
  lowPred = zeroes(lowPred)
  upPred = apply(preds, 2, up)

# Add a polygon for 95% CI
  polygon(c(Age, rev(Age)), c(lowPred, rev(upPred)),
          col=rgb(.8, .8, .8, 0.5),
          border=NA
          )

# Plot the mean and 95% CRI for predicted length at each age
  lines(Age, muPred, col='black', lwd=.9, lty=1)

# Plot axes
  axis(1, pos=0)
  axis(2, pos=0, las=2)

# Label axes
  mtext(expression(paste('Age (years)')), side=1, line=3)
  mtext('Total length (mm)', side=2, line=3.5, adj=1.2) 

# Close file connection
  dev.off()  


# . Figure 3. Covariate effects on Linf, K, and M -----
# .. Data ----
# Re-create the data used for Stan
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
tiff(filename='results/Figure3.tiff',
     width = 1300,
     height = 2000,
     res = 400,
     pointsize = 8
     )

# Set graphical parameters  
  par(mfrow=c(3,1), oma=c(4,5,1,1), mar=c(1,1,1,1))

# .. Hdrilla vs L-infinity -----
# Make a sequence of new biomasses
  newBiomass = seq(min(vb_data$hydrilla),
                   max(vb_data$hydrilla),
                   0.1
                   )

# Make an empty matrix to then fill with data later
  preds = matrix(NA, nrow = length(pars$b0_linf), ncol = length(newBiomass))

# Fill the matrix with predicted L-inf
# based on the linear predictor that
# uses hydrilla biomass
  for(i in 1:nrow(preds)){
    for(t in 1:length(newBiomass)){
      preds[i, t] = exp(pars$Gamma[i,1] + pars$b0_linf[i] + pars$bh_linf[i]*newBiomass[t])
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
  
# Add axes and labels  
  mtext(side=2, expression(paste(hat(italic('L'))[infinity])), line=4, cex=.66)
  axis(2, las = 2)
  axis(1, at=nscale(seq(0,1400,200), fish$ha), labels=FALSE)

# Calculate the mean and 95% CRIs for posterior predictions
  muPred = apply(preds, 2, mean)
  lowPred = apply(preds, 2, low)
  lowPred = zeroes(lowPred)
  upPred = apply(preds, 2, up)

# Plot 95% CRI polygon
  polygon(c(newBiomass, rev(newBiomass)),
          c(lowPred, rev(upPred)),
          col=rgb(.8, .8, .8, 0.5),
          border=NA
          )

# Plot posterior predictive mean
  lines(x = newBiomass, y = muPred, col = 'black', lwd = .5)


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
      preds[i, t] = exp(pars$Gamma[i,2] + pars$b0_k[i] + pars$bh_k[i]*newBiomass[t])
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
       ylim = c(0, 0.4)
       )
  
# Add axes and labels
  axis(1, at=nscale(seq(0,1400,200), fish$ha), labels=FALSE)
  mtext(side=2, expression(hat(italic('K'))), line=4, cex=.66)
  axis(2, las = 2)

# Calculate the mean and 95% CRIs for posterior predictions
  muPred = apply(preds, 2, mean)
  lowPred = apply(preds, 2, low)
  lowPred = zeroes(lowPred)
  upPred = apply(preds, 2, up)

# Plot 95% CRI polygon  
  polygon(c(newBiomass, rev(newBiomass)),
          c(lowPred, rev(upPred)),
          col=rgb(.8, .8, .8, 0.5),
          border=NA
          )

# Plot posterior predictive mean
  lines(x = newBiomass, y = muPred, col = 'black', lwd = .5)


# .. Hydrilla vs M -----
# Make an empty matrix to then fill with data later
  preds = matrix(NA, nrow = length(pars$b0_linf), ncol = length(newBiomass))
  
# Fill the matrix with predicted L-inf
# based on the linear predictor that
# uses hydrilla biomass
  for(i in 1:nrow(preds)){
    for(t in 1:length(newBiomass)){
      preds[i, t] = exp(pars$Gamma[i,2] + pars$b0_k[i] + pars$bh_k[i]*newBiomass[t])
    }
  }
  preds <- apply(preds, 2, FUN='*', 1.5)

# Plot the first prediction to get a plotting window set up
  plot(x = newBiomass, y = preds[1, ], type = 'l',
       col = rgb(0.6, 0.6, 0.6, 0.05),
       xlim = c(min(vb_data$hydrilla), max(vb_data$hydrilla)),
       xlab = '', xaxt='n',
       ylab = '',
       yaxt = 'n',
       ylim = c(0, .4)
       ) 
  
# Add axes and labels  
  mtext(side=1, "Surface hectares of hydrilla", line=3.5, cex=.66)
  mtext(side=2, expression(hat(italic('M'))), line=4, cex=.66)
  axis(1, at=nscale(seq(0,1400,200), fish$ha),
       labels=seq(0,1400,200))
  axis(2, at=seq(0.0,1,.1), labels=sprintf('%.1f', seq(0.0,1,.1)), las = 2)

# Calculate the mean and 95% CRIs for posterior predictions
  muPred = apply(preds, 2, mean)
  lowPred = apply(preds, 2, quantile, probs=0.025)
  upPred = apply(preds, 2, quantile, probs=0.975)

# Add 95% CRI polygon
polygon(c(newBiomass, rev(newBiomass)),
        c(lowPred, rev(upPred)),
        col=rgb(.8, .8, .8, 0.5),
        border=NA
        )

# Plot posterior predictive mean
  lines(x = newBiomass, y = muPred, col = 'black', lwd = .5)

# Close file connection
  dev.off()


# REVISIONS ----
# . Figure 1.R. Raw data plots ----
# Open file connection
  tiff(
    filename = paste0('results/Figure1.tiff'),
    height = 1500,
    width = 900,
    res = 400,
    pointsize = 6
  )
  
# Make the plot  
  par(mfrow=c(2,1), oma=c(4,4,0,0), mar=c(0,1,1,1))
  boxplot(Length ~ Year, data = fish,
          outline = FALSE, col='gray87',
          medlwd = 2, staplewex = 0, staplecol = NA,
          whisklty = 1, whisklwd = 2,
          boxlwd = 2, boxwex = 0.5,
          ylab = '', yaxt = 'n', ylim=c(0,1500), xaxt="n")
  axis(side = 1, tick = TRUE, labels = FALSE)
  axis(side = 2, las = 2)
  mtext('Total length (mm)', side=2, line = 3.5)
  stripchart(Length ~ Year, vertical = TRUE, data = fish, 
             method = "jitter", add = TRUE, pch = 20,
             col = 'black', cex=.75)
  boxplot(Age ~ Year, data = fish,
          outline = FALSE, col='gray87',
          medlwd = 2, staplewex = 0, staplecol = NA,
          whisklty = 1, whisklwd = 2,
          boxlwd = 2, boxwex = 0.5,
          ylab = '', yaxt = 'n',ylim=c(0,25))
  axis(side = 2, las = 2)
  mtext('Year of collection', side = 1, line=2.8)
  mtext('Age (years)', side=2, line = 3.5)
  stripchart(Age ~ Year, vertical = TRUE, data = fish, 
             method = "jitter", add = TRUE, pch = 20,
             col = 'black', cex=.75)
# Turn off graphics dev
  dev.off()
  

# . Figure 2.R. Growth curve ----
# Load the results
  load("results/vonbert_cohort-x-nha.rda")
  
# Extract parameter estimates from model fit
  pars = extract(fit)

# Create random draw from standard normal to
# include uncertainty related to standardized
# covariate
  # covs <- rnorm(length(pars$b0_k), 0, 1)
  covs <- sample(unique(scale(fish$nha)), length(pars$b0_k), replace = TRUE)

# covs <- runif(length(pars$b0_k), min(scale(fish$nha)), max(scale(fish$nha)))
  
# Offsets and covariate uncertainty
  Linf = exp(mean(pars$Gamma[,1,]) + pars$b0_linf + pars$ba_linf*covs)
  K = exp(mean(pars$Gamma[,2,]) + pars$b0_k + pars$ba_k*covs)
  t0 = exp(mean(pars$Gamma[,3,]) + pars$b0_t0)-10

# Make a sequence of new ages for which we will predict lengths
  Age = seq(1, 23, 1)

# Predict mean length at age for each sample
  preds = matrix(0, length(Linf), length(Age))
  for(i in 1:length(Linf)){
    for(t in 1:length(Age)){
      preds[i,t] = Linf[i]*(1-exp(-K[i]*(Age[t]-(t0[i]))))
    }
  }

# Set up graphics device
  tiff(
    filename = paste0('results/Figure2.tiff'),
    height = 1200,
    width = 1500,
    res = 400,
    pointsize = 6
  )
  
# Make the plot
  par(mfrow = c(1,1), mar = c(5,5,1,1))
  plot(fish$Age, fish$Length, 
       ylab = '',
       ylim=c(0, 1400),
       yaxt='n', 
       xlab='',
       xlim=c(0,25),
       axes = FALSE,
       pch = 21,
       col=rgb(.4, .4, .4, 0.5),
       bg=rgb(.4, .4, .4, 0.5),
       main='')

# Calculate the mean and 95% CRIs for posterior predictions
  muPred = apply(preds, 2, mean)
  lowPred = apply(preds, 2, quantile, probs=0.025)
  upPred = apply(preds, 2, quantile, probs=0.975)

# Add a polygon for 95% CI
  polygon(c(Age, rev(Age)), c(lowPred, rev(upPred)),
          col=rgb(.8, .8, .8, 0.5),
          border=NA
          )

# Plot the mean predicted length at each age
  lines(Age, muPred, col='black', lwd=.9, lty=1)
  
# Length at age for nha = 3 and nha = 275 
  # Standardized covariate values
    sx21 = (21 - mean(fish$nha))/sd(fish$nha)
    sx127 = (127 - mean(fish$nha))/sd(fish$nha)
  
  # Offsets and covariate uncertainty
    Linf21 = exp(mean(pars$Gamma[,1,]) + pars$b0_linf + pars$ba_linf*sx21)
    K21 = exp(mean(pars$Gamma[,2,]) + pars$b0_k + pars$ba_k*sx21)
    t021 = exp(mean(pars$Gamma[,3,]) + pars$b0_t0)-10
    Linf127 = exp(mean(pars$Gamma[,1,]) + pars$b0_linf + pars$ba_linf*sx127)
    K127 = exp(mean(pars$Gamma[,2,]) + pars$b0_k + pars$ba_k*sx127)
    t0127 = exp(mean(pars$Gamma[,3,]) + pars$b0_t0)-10  
    
  # Make a sequence of new ages for which we will predict lengths
    Age = seq(1, 23, 1)
    
    L21 = mean(Linf21)*(1-exp(-mean(K21)*(Age-(mean(t021)))))
    L127 = mean(Linf127)*(1-exp(-mean(K127)*(Age-(mean(t0127)))))

    lines(Age, L21, lty = 2)
    lines(Age, L127, lty = 2)    
    
# Plot axes
  axis(1, pos=0)
  axis(2, pos=0, las=2)

# Label axes
  mtext(expression(paste('Age (years)')), side=1, line=3)
  mtext('Total length (mm)', side=2, line=3.5)#, adj=1.2)
 
# Turn off graphics dev
  dev.off()
  
# . Figure 3.R. Covariate effects on Linf, K, and M ----- 
# .. Data ----
fish = fish[fish$agec == fish$Age, ]  
  
# .. Load the results ----
load("results/vonbert_cohort-x-nha.rda")
  
# Extract posteriors from model object
  pars = extract(fit)

# .. Set up an image file ----
tiff(filename='results/Figure3.tiff',
     width = 1300,
     height = 2000,
     res = 400,
     pointsize = 8
     )

# Set graphical parameters  
  par(mfrow=c(3,1), oma=c(4,5,1,1), mar=c(1,1,1,1))
  
# .. nha vs Linf ----  
# Make a new vector of scaled abundances
  n = seq(min(scale(fish$nha))-.05, max(scale(fish$nha))+.3, 0.3)

# Make scaled abundances for plotting
  new = seq(0, round(max(fish$nha)+3e2, 1), 25)
  new
  scaled.new <- scale(new, center = mean(fish$nha), scale = sd(fish$nha))
    
# Predict Linf for each sample
  lpreds = matrix(data = NA, nrow=length(pars$b0_linf), ncol=length(n))
  for(i in 1:length(pars$b0_linf)){
    for(t in 1:length(n)){
      lpreds[i, t] = exp(pars$b0_linf[i] + pars$ba_linf[i]*n[t])
    }
  }
  
# Calculate the mean and 95% CRIs for posterior predictions
  muPred = apply(lpreds, 2, mean)
  lowPred = apply(lpreds, 2, quantile, probs=0.025)
  upPred = apply(lpreds, 2, quantile, probs=0.975)  

# Summarize individual posteriors for Linf for plotting
  linfs <- data.frame(
    linfmeans = apply(pars$Linf, 2, mean),
    linfl = apply(pars$Linf, 2, quantile, probs=0.025),
    linfu = apply(pars$Linf, 2, quantile, probs=0.975)
  )
  
# Jitter the scaled covariate values and keep only 
# the unique combos of cohort x nha (unique Linfs)
  ns <- scale(fish$nha)#[!duplicated(linfs[,1])]
  ns <- jitter(ns, 0)
  #linfs <- linfs[!duplicated(linfs[,1]), ]

# Add points
  plot(x = ns,
       y = linfs[,1],
       pch=21, bg='black',
       ylim=c(1000, 1500),
       yaxt='n', xlab='', ylab='',
       xlim=c(min(n), round(max(n), 1)), axes = F,
       main='')   
  

# Add a polygon for 95% CI
  polygon(c(n, rev(n)), c(lowPred, rev(upPred)),
          col=rgb(.8, .8, .8, 0.25),
          border=rgb(.8, .8, .8, 0.25)
  )

# Add 95% CRI as jittered line segments
  segments(x0=ns,
           y0=linfs[,2],
           y1=linfs[,3],
           col=rgb(0,0,0,0.05), lwd=.5)
  points(x=ns, y=linfs[,1], pch=21,
         bg=rgb(0,0,0,0.01),
         col=rgb(0,0,0,0.01)
         ) #grey.colors(5)[as.factor(fish$Year)])
           #'black', col='black')
  
# Plot posterior predictive mean
  lines(n, muPred, col='gray40', lwd=1, lty=1)
  
# Add axes and labels
  axis(1, pos=1000, at=scaled.new, tick=TRUE, labels=FALSE)
  axis(2, pos=min(n)-0.05, las=2)  
  mtext(expression(
    paste(italic('L')[infinity])),
    side=2, line=3.5
    )
  
 
# .. nha vs K ----
  # Predict K for each sample
  kpreds = matrix(data = NA, nrow=length(pars$b0_k), ncol=length(n))
  for(i in 1:length(pars$b0_k)){
    for(t in 1:length(n)){
      kpreds[i, t] = exp(pars$b0_k[i] + pars$ba_k[i]*n[t])
    }
  }
  
# Calculate the mean and 95% CRIs for posterior predictions
  muPred = apply(kpreds, 2, mean)
  lowPred = apply(kpreds, 2, quantile, probs=0.025)
  upPred = apply(kpreds, 2, quantile, probs=0.975)  

# Summarize individual posteriors for k for plotting
  ks <- data.frame(
    kmeans = apply(pars$K, 2, mean),
    kl = apply(pars$K, 2, quantile, probs=0.025),
    ku = apply(pars$K, 2, quantile, probs=0.975)
  )
  
# Jitter the scaled covariate values and keep only 
# the unique combos of cohort x nha (unique ks)
  ns <- scale(fish$nha)#[!duplicated(ks[,1])]
  ns <- jitter(ns, 0)
  #ks <- ks[!duplicated(ks[,1]), ]

# Add points
  plot(x = ns,
       y = ks[,1],
       pch=21, bg='black',
       ylim=c(0, .3),
       yaxt='n', xlab='', ylab='',
       xlim=c(min(n), round(max(n), 1)), axes = F,
       main='')   
  

# Add a polygon for 95% CI
  polygon(c(n, rev(n)), c(lowPred, rev(upPred)),
          col=rgb(.8, .8, .8, 0.25),
          border=rgb(.8, .8, .8, 0.25)
  )

# Add 95% CRI as jittered line segments
  segments(x0=ns,
           y0=ks[,2],
           y1=ks[,3],
           col=rgb(0,0,0,0.05), lwd=.5)
  points(x=ns, y=ks[,1], pch=21,
         bg=rgb(0,0,0,0.01),
         col=rgb(0,0,0,0.01))
  
# Plot posterior predictive mean
  lines(n, muPred, col='gray40', lwd=1, lty=1)
  
# Add axes and labels
  axis(1, pos=0, at=scaled.new, tick=TRUE, labels=FALSE)
  axis(2, pos=min(n)-0.05, las=2)  
  mtext(expression(
    paste(italic('K'))),
    side=2, line=3.5
    )
  
# .. nha vs M ----
# Predict M for each sample as M = 4.118 * (K^0.73) * (L^−0.33) 
  preds = matrix(data = NA, nrow=length(pars$b0_k), ncol=length(n))
  for(i in 1:length(pars$b0_k)){
    for(t in 1:length(n)){
      preds[i, t] = 4.118 * (kpreds[i,t]^0.73) * ((lpreds[i,t]/10)^-0.33)
    }
  }
  
# Calculate the mean and 95% CRIs for posterior predictions
  muPred = apply(preds, 2, mean)
  lowPred = apply(preds, 2, quantile, probs=0.025)
  upPred = apply(preds, 2, quantile, probs=0.975)  

# Summarize individual posteriors for M for plotting
  ms <- data.frame(
    mmeans = 4.118*(ks[,1]^0.73)*(linfs[,1]/10)^(-0.33),
    ml = 4.118*(ks[,2]^0.73)*(linfs[,2]/10)^(-0.33),
    mu = 4.118*(ks[,3]^0.73)*(linfs[,3]/10)^(-0.33)
  )
  
# Jitter the scaled covariate values and keep only 
# the unique combos of cohort x nha (unique ks)
  ns <- scale(fish$nha)#[!duplicated(ks[,1])]
  ns <- jitter(ns, 0)
  #ks <- ks[!duplicated(ks[,1]), ]

# Add points
  plot(x = ns,
       y = ms[,1],
       pch=21, bg='black',
       ylim=c(0, .5),
       yaxt='n', xlab='', ylab='',
       xlim=c(min(n), round(max(n), 1)), axes = F,
       main='')
  

# Add a polygon for 95% CI
  polygon(c(n, rev(n)), c(lowPred, rev(upPred)),
          col=rgb(.8, .8, .8, 0.25),
          border=rgb(.8, .8, .8, 0.25)
  )

# Add 95% CRI as jittered line segments
  segments(x0=ns,
           y0=ms[,2],
           y1=ms[,3],
           col=rgb(0,0,0,0.05), lwd=.5)
  points(x=ns, y=ms[,1], pch=21,
         bg=rgb(0,0,0,0.01),
         col=rgb(0,0,0,0.01))
  
# Plot posterior predictive mean
  lines(n, muPred, col='gray40', lwd=1, lty=1)
  
# Add axes and labels
  axis(1, pos=0, at=scaled.new, labels=new)
  axis(2, pos=min(n)-0.05, las=2)  
  mtext(expression(
    paste('Grass Carp density (n ', '\u00b7', ' ha'^'-1',')')),
    side=1, line=3.5
    )
  mtext(expression(
    paste(italic('M'))),
    side=2, line=3.5
    )
  
# Turn off graphics device
dev.off()  
  
  