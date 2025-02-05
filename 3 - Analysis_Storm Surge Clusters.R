##########################################################################
#
# Analysis of clusters of storm closures
# Copyright (C) 2025  Alexander Bakker (TU Delft/Rijkswaterstaat)
#
# This script is developed as part of the study
#  "Storm Surge Clusters, Multi-Peaked Storms and Their Effect on
#   the Performance of the Maeslant Storm Surge Barrier (the  Netherlands)." 
#   Bakker, A.M.R.; Rovers, D.L.T.; Mooyaart, L.F. (2025)
#   J. Mar. Sci. Eng. 2025,13, x. https://doi.org/10.3390/xxxxx
#
# In this script the peak water levels that exceed 
# MSL +212.2cm (MSL +xxx.xcm after bias correction) are analysed
# as described in subsections 3.1 and 4.1
# 
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
###########################################################################

###################################################
#   Initializing                                  #
###################################################

rm(list=ls())

require(extRemes)
require(evmix)

# input files
file="Data/SLRpeaks.csv"

# input variables
dpy   = 365.25     # number of days per year
hpd   = 24         # number of hours per day
thr   = 212.2      # applied threshold [cmMSL]
nyr   = 65         # number of years (1953-2018)
cdl   = 300        # closure decision level [cmMSL]
cfr   = 0.1        # (official) closure frequency [yr-1]

# experiments
slrs  = data.frame(exp=c('slr0','slr25','slr50'),  # sea level rise
                   slr=c(0,25,50))

# Read data
data        <- read.table(file, header=TRUE, sep=";", dec=".")
names(data) <- c("dt","SL")
data$dt     <- as.POSIXct(data$dt, format = "%d-%m-%Y %H:%M")
n           <- length(data$dt) # length of data set
fr          <- n/nyr           # number of extreme surges exceeding threshold

###################################################
#   Analysis interarrival times                   #
###################################################
a  <- as.numeric((data$dt[2:n] - data$dt[1:n-1])   ) # interarrivel times

# randomly generated interarrival times accounting for seasonality
frm = rep(                                           # peak prob per day
  c(rep(25/65/31,31),
    rep(27/65/28,28),
    rep(10/65/31,31),
    rep(3/65/30,30),
    rep(0,31),
    rep(0,30),
    rep(0,31),
    rep(0,31),
    rep(4/65/30,30),
    rep(14/65/31,31),
    rep(30/65/30,30),
    rep(30/65/31,31)),
  10000)
sd <- which(frm > runif(length(frm)))              # randomly generated storm days (accounting for seasonality)
rn <- length(sd)
ra <- sd[2:rn] - sd[1:rn-1]        # randomly generated interarrival times

# estimate CDF interarrival times
nP     <- 1050               # number of days analysed
cdfa   <- rep(NA,nP)
ecdfa  <- rep(NA,nP)
recdfa <- rep(NA,nP)
lambda <- n/(365.25 * nyr)   # lambda theoretical CDF

for (i in 1:nP) {
  ecdfa[i]  <- sum(a<=i)/(n-1)    # empirical cdf interarrival times observations
  recdfa[i] <- sum(ra<=i)/(rn-1)  # empirical cdf interarrival times generated data
  cdfa[i]   <- 1-exp(-lambda*i)   # theoretical cdf interarrival times
}

# plot CDF interarrival times
matplot(cdfa, type='l', col='black',lwd=3,
        cex.lab=1.4,cex.axis=1.2,
        xlab='Number of days',
        ylab='CDF')#,
lines(recdfa, col='gray',lwd=3)
lines(ecdfa, col='red',lwd=3)
legend("bottomright", inset = 0.05,
       legend = c("Theoretical",
                  "Theoretical (with seasonality)",
                  "Empirical"),
       bg  = 'white', # Legend background color
       box.col = "white",
       lty = c(1,1,1),
       col = c('black','gray','red'),
       lwd = 3,
       cex = 1.4)

# compare
head(table(ceiling(a))/65,25)
head(table(ceiling(ra))/10000,25)

rat <- as.data.frame(head(table(ceiling(ra))/10000,1000))
plot(rat$Var1,rat$Freq)

# check dependence peaks and interarrival time
x1 <- data$SL[1:n-1]  # first peak in sequence of two
x2 <- data$SL[2:n]    # second peak in sequence of two
dx <- x1 - x2          # difference of two successive peaks
mx <- (x1+x2)/2        # mean of two successive peaks

cor(a,x1,method="spearman")  # correlation interarrival time and first peak
cor(a,x2,method="spearman")  # correlation interarrival time and second peak
cor(a,dx,method="spearman")  # correlation interarrival time and difference peaks
cor(a,mx,method="spearman")  # correlation interarrival time and mean of peaks
cor(x1,x2,method="spearman") # correlation first and second peak

# compare stats of peaks of twin storm within two days
x1a2 <- x1[which(a<=2)]         # select first peaks
x2a2 <- x2[which(a<=2)]         # select second peaks

c(mean(x1),   sd(x1),   quantile(x1,  c(0.25,0.5,0.75)))
c(mean(x2),   sd(x2),   quantile(x2,  c(0.25,0.5,0.75)))
c(mean(x1a2), sd(x1a2), quantile(x1a2,c(0.25,0.5,0.75)))
c(mean(x2a2), sd(x2a2), quantile(x2a2,c(0.25,0.5,0.75)))

#########################################################
#   Fit GPD to sea water level peaks and determine bias #
#########################################################
# fit GPD (! bias correction won't change scale and shape parameter)
fitD  <- fevd(SL, data, threshold = thr, type = "GP",
              time.units = paste0(fr,"/year"))

# extract GPD parameters
scale  <- fitD$results$par['scale']
shape  <- fitD$results$par['shape']

# determine bias
bias <- cdl - qgpd(1-(cfr/fr), u=thr, sigmau=scale, xi=shape, phiu=1)

###################################################
#   Conditional exceedance probabilities &        #
#           interarrival times closure events     #
###################################################
# output/variables
freq      <- data.frame()        # closure frequency
vals      <- seq(300,700,1)
cp.single <- data.frame(hx=vals) # conditional probabilities

# estimate emperical discrete pdf interarrival time storm events
epdfa     <- rep(NA,nP)
for (i in 1:nP) {
  epdfa[i]  <- sum(ceiling(a)==i)/(n-1)   # empirical pdf interarrival times observations
}

# sample interarrival times closure events
nsamples  <- 100000                             # number of storm surge samples
tsamples  <- cumsum(sample.int(n=nP,            # timing storm events
                               size=nsamples,
                               replace=T,
                               prob=epdfa))

# start experiments
for (i in 1:nrow(slrs)) {
  threshold <- thr + bias + slrs$slr[i]         # bias corrected threshold for particular slr
  hxsamples <- rgpd(n=nsamples, u=threshold,sigmau=scale,xi=shape,phiu=1) # random sampling
  
  # conditional closure probability (given a storm event) and resulting Closure frequency
  Pcdl     <- (1 - pgpd(cdl,u=threshold,sigmau=scale,xi=shape,phiu=1)) # conditional closure probability
  Fcl      <- fr * Pcdl                                                # closure frequency
  
  # interarrival times closure events
  asamples <- tsamples[which(hxsamples>300)]
  asamples <- asamples[2:length(asamples)] - asamples[1:length(asamples)-1]
  
  # probability that previous closure was within certain time frame before
  Pcl.2d   <-                      sum(asamples<=2)  / length(asamples)
  Pcl.w    <- (sum(asamples<=7)  - sum(asamples<=2)) / length(asamples)
  Pcl.m    <- (sum(asamples<=30) - sum(asamples<=7)) / length(asamples)
  Pcl.mp   <- 1 - Pcl.2d - Pcl.w - Pcl.m
  
  freq  <- rbind(freq, data.frame(
             exp = slrs$exp[i],
             slr = slrs$slr[i],
              fr = Fcl,
        two.days = Pcl.2d,
            week = Pcl.w,
           month = Pcl.m,
      month.plus = Pcl.mp))

  # estimate conditional probabilities  
  ccdf     <- (pgpd(vals, u=threshold, sigmau=scale, xi=shape, phiu=1) -
               pgpd(cdl,  u=threshold, sigmau=scale, xi=shape, phiu=1) ) /
          (1 - pgpd(cdl,  u=threshold, sigmau=scale, xi=shape, phiu=1) )
  cp.single[slrs$exp[i]] <- 1 -  ccdf
}

write.csv(freq,"Results/frequency_closure_events_ssc_analysis.csv",row.names=F)
write.csv(cp.single,"Results/cond_exc_prob_single_ssc_analysis.csv",row.names=F)


# return periods
RP   <- pgpd(c(thr,213:400)+bias,    u = thr+bias,    sigmau = scale, xi = shape, phiu = 1)
eRP  <- rank(data$SL+bias)/(length(data$SL) + 1)

RP   <- 1 / (1 - RP)   / fr
eRP  <- 1 / (1 - eRP) / fr

plot(RP,c(thr,213:400)+bias,
     type='l', lwd=4, col='red',
     xlab="Return Period", ylab="Sea water level [cmMSL]",
     xaxt="n",
     log='x',
     xlim=c(0.1,1000),
     ylim=c(200,500),
     cex.lab=1.4,cex.axis=1.2)

lines(RP,c(thr,213:400)+bias+25,
      type='l',
      lwd=4,
      col='lightblue')

lines(RP,c(thr,213:400)+bias+50,
      type='l',
      lwd=4,
      col='blue')


lines(eRP,data$SL+bias,
      type='p',
      pch=16,
      cex=1.2,
      col='black')

lines(eRP,data$SL,
      type='p',
      pch=16,
      cex=1.2,
      col='gray')


abline(h=c(300,360), lty=2, lwd=1.5 )
text(x=0.07, y=310, pos=4, cex=1.2, "closure decision level (MSL +300cm)")
text(x=0.07, y=370, pos=4, cex=1.2, "critical level interior flood defences (MSL +360cm)")

axis(1, at = c(0.1, 1, 10, 100, 1000, 10000),
     labels = c("0.1", "1", "10", "100", "1000", "10000"),
     cex.axis=1.2)
legend("topleft", inset = 0.01,
       legend = c(' 0 cm sea level rise',
                  '25 cm sea level rise',
                  '50 cm sea level rise',
                  'observations (corrected)',
                  'observations (uncorrected)'),
       bg = 'white', # Legend background color
       box.col = "white",
       pch = c(NA,NA,NA,16,16),
       lty = c(1,1,1,NA,NA),
       col = c('red','lightblue','blue','black','gray'),
       lwd = c(2,2,2,0,0),
       cex = 1.2)
 



