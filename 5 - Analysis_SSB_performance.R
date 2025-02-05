##########################################################################
#
# Analysis of Storm Surge Barrier Performance
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

### load packages
#require(extRemes)
#require(evmix)
#require(plotrix)

### input variables
# Probability of non closure
Pf.ic = 0.01                    # base/isolated closure
Pf.pc = c( 0.01,  0.05,   0.1)  # preceded closure
Pf.sc = c( 0.001, 0.005, 0.01)  # second closure

### read files
# frequency tables
ssc.freq <- read.csv('Results/frequency_closure_events_ssc_analysis.csv', header=T, sep=',')
mps.freq <- read.csv('Results/frequency_closure_events_mps_analysis.csv', header=T, sep=';')

# conditional exceedance probabilities
#cp.single <- read.csv('Results/cond_exc_prob_single_ssc_analysis.csv', header=T, sep=',')
cp.single <- read.csv('Results/cond_exc_prob_single_mps_analysis.csv', header=T, sep=',')
cp.double <- read.csv('Results/cond_exc_prob_double.csv', header=T, sep=',')
cp.first  <- read.csv('Results/cond_exc_prob_first.csv', header=T, sep=',')
cp.second <- read.csv('Results/cond_exc_prob_second.csv', header=T, sep=',')

####################
### experiment 1 ###
slr = 0
hin = 100

# probability that closure event is double closure
P.double   <- mps.freq$proportion_double[1]

# Effective probability of failed first or single closure
P.pc    <- 0.023                  # probability of preceding closure
Pf.fc   <- (1 - P.pc) * Pf.ic +   # effective failure probability
                P.pc  * Pf.pc

# closure.frequencies
cf.total   <- ssc.freq[1,3]
cf.single  <- cf.total * (1 - P.double)
cf.double  <- cf.total * (P.double)

# exceedance frequencies
hx             <- 300:700

#F.nobarrier     <- cf.total  * cp.single$slr0 # without barrier
#F.single.nb     <- cf.single * cp.single$slr0
F.single.nb     <- cf.single * cp.single$exp1 # without barrier
F.first.nb      <- cf.double * cp.first$exp1
F.second.nb     <- cf.double * cp.second$exp1
F.double.nb     <- cf.double * cp.double$exp1

F.nobarrier     <- F.single.nb + F.double.nb

F.default       <- Pf.ic * F.nobarrier        # default reliability
F.single.def    <- Pf.ic * F.single.nb
F.first.def     <- Pf.ic * F.first.nb
F.second.def    <- Pf.ic * F.second.nb

F.cluster       <- Pf.fc[2] * F.nobarrier     # effect of clustering
F.cluster.ll    <- Pf.fc[1] * F.nobarrier
F.cluster.ul    <- Pf.fc[3] * F.nobarrier

F.second.mps    <- Pf.sc[2] * F.second.nb     # effect of mps
F.second.mps.ll <- Pf.sc[1] * F.second.nb
F.second.mps.ul <- Pf.sc[3] * F.second.nb

F.mps           <- F.single.def + F.first.def + F.second.mps
F.mps.ll        <- F.single.def + F.first.def + F.second.mps.ll
F.mps.ul        <- F.single.def + F.first.def + F.second.mps.ul

F.single.ssc    <- Pf.fc[2] * F.single.nb     # combined effect
F.single.ssc.ll <- Pf.fc[1] * F.single.nb
F.single.ssc.ul <- Pf.fc[3] * F.single.nb
F.first.ssc     <- Pf.fc[2] * F.first.nb
F.first.ssc.ll  <- Pf.fc[1] * F.first.nb
F.first.ssc.ul  <- Pf.fc[3] * F.first.nb

F.combined      <- F.single.ssc    + F.first.ssc    + F.second.mps
F.combined.ll   <- F.single.ssc.ll + F.first.ssc.ll + F.second.mps.ll
F.combined.ul   <- F.single.ssc.ul + F.first.ssc.ul + F.second.mps.ul

# return periods
RP.nobarrier    <- 1 / F.nobarrier
RP.default      <- 1 / F.default
RP.cluster      <- 1 / F.cluster
RP.cluster.ll   <- 1 / F.cluster.ll
RP.cluster.ul   <- 1 / F.cluster.ul
RP.mps          <- 1 / F.mps
RP.mps.ll       <- 1 / F.mps.ll
RP.mps.ul       <- 1 / F.mps.ul
RP.combined     <- 1 / F.combined
RP.combined.ll  <- 1 / F.combined.ll
RP.combined.ul  <- 1 / F.combined.ul


# plot return periods
charttitle = "No sea level rise, initial inner water level = 100cm"

plot(RP.nobarrier,hx,
     type='l', lwd=3, lty=1, col='grey',
     xlab="Return Period [years]", ylab="Sea water level [cmMSL]",
     xaxt="n",
     log='x',
     xlim=c(1,1000000),
     ylim=c(300,550),
     main=charttitle,
     cex.lab=1.4,cex.axis=1.2)

abline(h=c(300,360), lty=2, lwd=1.5 )
text(x=1500000, y=307, pos=2, cex=1.1, "(MSL +300cm)")
text(x=1500000, y=367, pos=2, cex=1.1, "(MSL +360cm)")
axis(1, at = c(10, 100, 1000, 10000,100000),
     labels = c("10", "100", "1000", "10,000","100,000"),
     cex.axis=1.2)
axis(1, at = c(1:10,10*(1:10), 100*(1:10), 1000*(1:10),10000*(1:10),100000*(1:10)),
     labels=F,lwd.ticks=0.5,
     cex.axis=1.2)
legend("topleft", inset = 0.002,
       legend = c('without storm surge barrier',
                  'constant probability of non-closure',
#                  'effect on storm clustering',
#                  'effect of multi-closure events',
                  'combined effect of clustering and multi-closure events',
                  'upper and lower limits'),
       bg = 'white', # Legend background color
       box.col = "white",
       lty = c(1,1,1,2), #c(1,1,1,1,1,2),
#       col = c('grey','black','green3','turquoise','deepskyblue3','black'),
       col = c('grey','black','deepskyblue3','deepskyblue3'),
       lwd = c(3,3,3,1), #c(3,3,3,3,3,1),
       cex = 1.1)

lines(RP.default,hx, # default reliability
      type='l',
      lwd=3,
      col='black')

lines(RP.cluster,hx, # effect clustering
      type='l',
      lwd=3,
      col='green3')

lines(RP.cluster.ll,hx,
      type='l',
      lty=2,
      lwd=1,
      col='green3')

lines(RP.cluster.ul,hx,
      type='l',
      lty=2,
      lwd=1,
      col='green3')

lines(RP.mps,hx, # effect mps
      type='l',
      lwd=3,
      col='turquoise')

lines(RP.mps.ll,hx,
      type='l',
      lty=2,
      lwd=1,
      col='turquoise')

lines(RP.mps.ul,hx,
      type='l',
      lty=2,
      lwd=1,
      col='turquoise')

lines(RP.combined,hx, # combined effect
      type='l',
      lwd=3,
      col='deepskyblue3')

lines(RP.combined.ll,hx,
      type='l',
      lty=2,
      lwd=1,
      col='deepskyblue3')

lines(RP.combined.ul,hx,
      type='l',
      lty=2,
      lwd=1,
      col='deepskyblue3')

### end experiment 1 ###
########################


####################
### experiment 9 ###
slr = 50
hin = 100

# probability that closure event is double closure
P.double   <- mps.freq$proportion_double[9]

# Effective probability of failed first or single closure
P.pc    <- 0.184                  # probability of preceding closure
Pf.fc   <- (1 - P.pc) * Pf.ic +   # effective failure probability
                P.pc  * Pf.pc

# closure.frequencies
cf.total   <- ssc.freq[3,3]
cf.single  <- cf.total * (1 - P.double)
cf.double  <- cf.total * (P.double)

# exceedance frequencies
hx             <- 300:700

#F.nobarrier     <- cf.total  * cp.single$slr0 # without barrier
#F.single.nb     <- cf.single * cp.single$slr0
F.single.nb     <- cf.single * cp.single$exp9 # without barrier
F.first.nb      <- cf.double * cp.first$exp9
F.second.nb     <- cf.double * cp.second$exp9
F.double.nb     <- cf.double * cp.double$exp9

F.nobarrier     <- F.single.nb + F.double.nb

F.default       <- Pf.ic * F.nobarrier        # default reliability
F.single.def    <- Pf.ic * F.single.nb
F.first.def     <- Pf.ic * F.first.nb
F.second.def    <- Pf.ic * F.second.nb

F.cluster       <- Pf.fc[2] * F.nobarrier     # effect of clustering
F.cluster.ll    <- Pf.fc[1] * F.nobarrier
F.cluster.ul    <- Pf.fc[3] * F.nobarrier

F.second.mps    <- Pf.sc[2] * F.second.nb     # effect of mps
F.second.mps.ll <- Pf.sc[1] * F.second.nb
F.second.mps.ul <- Pf.sc[3] * F.second.nb

F.mps           <- F.single.def + F.first.def + F.second.mps
F.mps.ll        <- F.single.def + F.first.def + F.second.mps.ll
F.mps.ul        <- F.single.def + F.first.def + F.second.mps.ul

F.single.ssc    <- Pf.fc[2] * F.single.nb     # combined effect
F.single.ssc.ll <- Pf.fc[1] * F.single.nb
F.single.ssc.ul <- Pf.fc[3] * F.single.nb
F.first.ssc     <- Pf.fc[2] * F.first.nb
F.first.ssc.ll  <- Pf.fc[1] * F.first.nb
F.first.ssc.ul  <- Pf.fc[3] * F.first.nb

F.combined      <- F.single.ssc    + F.first.ssc    + F.second.mps
F.combined.ll   <- F.single.ssc.ll + F.first.ssc.ll + F.second.mps.ll
F.combined.ul   <- F.single.ssc.ul + F.first.ssc.ul + F.second.mps.ul

# return periods
RP.nobarrier    <- 1 / F.nobarrier
RP.default      <- 1 / F.default
RP.cluster      <- 1 / F.cluster
RP.cluster.ll   <- 1 / F.cluster.ll
RP.cluster.ul   <- 1 / F.cluster.ul
RP.mps          <- 1 / F.mps
RP.mps.ll       <- 1 / F.mps.ll
RP.mps.ul       <- 1 / F.mps.ul
RP.combined     <- 1 / F.combined
RP.combined.ll  <- 1 / F.combined.ll
RP.combined.ul  <- 1 / F.combined.ul


# plot return periods
charttitle = "50cm sea level rise, initial inner water level = 100cm"

plot(RP.nobarrier,hx,
     type='l', lwd=3, lty=1, col='grey',
     xlab="Return Period [years]", ylab="Sea water level [cmMSL]",
     xaxt="n",
     log='x',
     xlim=c(1,1000000),
     ylim=c(300,550),
     main=charttitle,
     cex.lab=1.4,cex.axis=1.2)

abline(h=c(300,360), lty=2, lwd=1.5 )
text(x=1500000, y=307, pos=2, cex=1.1, "(MSL +300cm)")
text(x=1500000, y=367, pos=2, cex=1.1, "(MSL +360cm)")
axis(1, at = c(10, 100, 1000, 10000,100000),
     labels = c("10", "100", "1000", "10,000","100,000"),
     cex.axis=1.2)
axis(1, at = c(1:10,10*(1:10), 100*(1:10), 1000*(1:10),10000*(1:10),100000*(1:10)),
     labels=F,lwd.ticks=0.5,
     cex.axis=1.2)
legend("topleft", inset = 0.002,
       legend = c('without storm surge barrier',
                  'constant probability of non-closure',
                  #'effect on storm clustering',
                  #'effect of multi-closure events',
                  'combined effect of clustering and multi-closure events',
                  'upper and lower limits'),
       bg = 'white', # Legend background color
       box.col = "white",
       lty = c(1,1,1,2), #c(1,1,1,1,1,2),
       #col = c('grey','black','green3','turquoise','deepskyblue3','black'),
       col = c('grey','black','deepskyblue3','deepskyblue3'),
       lwd = c(3,3,3,1), #c(3,3,3,3,3,1),
       cex = 1.1)

lines(RP.default,hx, # default reliability
      type='l',
      lwd=3,
      col='black')

lines(RP.cluster,hx, # effect clustering
      type='l',
      lwd=3,
      col='green3')

lines(RP.cluster.ll,hx,
      type='l',
      lty=2,
      lwd=1,
      col='green3')

lines(RP.cluster.ul,hx,
      type='l',
      lty=2,
      lwd=1,
      col='green3')

lines(RP.mps,hx, # effect mps
      type='l',
      lwd=3,
      col='turquoise')

lines(RP.mps.ll,hx,
      type='l',
      lty=2,
      lwd=1,
      col='turquoise')

lines(RP.mps.ul,hx,
      type='l',
      lty=2,
      lwd=1,
      col='turquoise')

lines(RP.combined,hx, # combined effect
      type='l',
      lwd=3,
      col='deepskyblue3')

lines(RP.combined.ll,hx,
      type='l',
      lty=2,
      lwd=1,
      col='deepskyblue3')

lines(RP.combined.ul,hx,
      type='l',
      lty=2,
      lwd=1,
      col='deepskyblue3')

### end experiment 9 ###
########################


####################
### experiment 12+ ###
slr = 50
hin = 100

Pf.sc = c( 0.01, 0.05, 0.1)

# probability that closure event is double closure
P.double   <- mps.freq$proportion_double[12]

# Effective probability of failed first or single closure
P.pc    <- 0.184                  # probability of preceding closure
Pf.fc   <- (1 - P.pc) * Pf.ic +   # effective failure probability
                P.pc  * Pf.pc

# closure.frequencies
cf.total   <- ssc.freq[3,3]
cf.single  <- cf.total * (1 - P.double)
cf.double  <- cf.total * (P.double)

# exceedance frequencies
hx             <- 300:700

#F.nobarrier     <- cf.total  * cp.single$slr0 # without barrier
#F.single.nb     <- cf.single * cp.single$slr0
F.single.nb     <- cf.single * cp.single$exp12 # without barrier
F.first.nb      <- cf.double * cp.first$exp12
F.second.nb     <- cf.double * cp.second$exp12
F.double.nb     <- cf.double * cp.double$exp12

F.nobarrier     <- F.single.nb + F.double.nb

F.default       <- Pf.ic * F.nobarrier        # default reliability
F.single.def    <- Pf.ic * F.single.nb
F.first.def     <- Pf.ic * F.first.nb
F.second.def    <- Pf.ic * F.second.nb

F.cluster       <- Pf.fc[2] * F.nobarrier     # effect of clustering
F.cluster.ll    <- Pf.fc[1] * F.nobarrier
F.cluster.ul    <- Pf.fc[3] * F.nobarrier

F.second.mps    <- Pf.sc[2] * F.second.nb     # effect of mps
F.second.mps.ll <- Pf.sc[1] * F.second.nb
F.second.mps.ul <- Pf.sc[3] * F.second.nb

F.mps           <- F.single.def + F.first.def + F.second.mps
F.mps.ll        <- F.single.def + F.first.def + F.second.mps.ll
F.mps.ul        <- F.single.def + F.first.def + F.second.mps.ul

F.single.ssc    <- Pf.fc[2] * F.single.nb     # combined effect
F.single.ssc.ll <- Pf.fc[1] * F.single.nb
F.single.ssc.ul <- Pf.fc[3] * F.single.nb
F.first.ssc     <- Pf.fc[2] * F.first.nb
F.first.ssc.ll  <- Pf.fc[1] * F.first.nb
F.first.ssc.ul  <- Pf.fc[3] * F.first.nb

F.combined      <- F.single.ssc    + F.first.ssc    + F.second.mps
F.combined.ll   <- F.single.ssc.ll + F.first.ssc.ll + F.second.mps.ll
F.combined.ul   <- F.single.ssc.ul + F.first.ssc.ul + F.second.mps.ul

# return periods
RP.nobarrier    <- 1 / F.nobarrier
RP.default      <- 1 / F.default
RP.cluster      <- 1 / F.cluster
RP.cluster.ll   <- 1 / F.cluster.ll
RP.cluster.ul   <- 1 / F.cluster.ul
RP.mps          <- 1 / F.mps
RP.mps.ll       <- 1 / F.mps.ll
RP.mps.ul       <- 1 / F.mps.ul
RP.combined     <- 1 / F.combined
RP.combined.ll  <- 1 / F.combined.ll
RP.combined.ul  <- 1 / F.combined.ul


# plot return periods
charttitle = "50cm sea level rise, initial inner water level = 100cm"

plot(RP.nobarrier,hx,
     type='l', lwd=3, lty=1, col='grey',
     xlab="Return Period [years]", ylab="Sea water level [cmMSL]",
     xaxt="n",
     log='x',
     xlim=c(1,1000000),
     ylim=c(300,550),
     main=charttitle,
     cex.lab=1.4,cex.axis=1.2)

abline(h=c(300,360), lty=2, lwd=1.5 )
text(x=1500000, y=307, pos=2, cex=1.1, "(MSL +300cm)")
text(x=1500000, y=367, pos=2, cex=1.1, "(MSL +360cm)")
axis(1, at = c(10, 100, 1000, 10000,100000),
     labels = c("10", "100", "1,000", "10,000","100,000"),
     cex.axis=1.2)
axis(1, at = c(1:10,10*(1:10), 100*(1:10), 1000*(1:10),10000*(1:10),100000*(1:10)),
     labels=F,lwd.ticks=0.5,
     cex.axis=1.2)
legend("topleft", inset = 0.002,
       legend = c('without storm surge barrier',
                  'constant probability of non-closure',
                  'effect on storm clustering',
                  'effect of multi-closure events',
                  'combined effect of clustering and multi-closure events',
                  'upper and lower limits'),
       bg = 'white', # Legend background color
       box.col = "white",
       lty = c(1,1,1,1,1,2),
       col = c('grey','black','green3','turquoise','deepskyblue3','black'),
       lwd = c(3,3,3,3,3,1),
       cex = 1.1)

lines(RP.default,hx, # default reliability
      type='l',
      lwd=3,
      col='black')

lines(RP.cluster,hx, # effect clustering
      type='l',
      lwd=3,
      col='green3')

lines(RP.cluster.ll,hx,
      type='l',
      lty=2,
      lwd=1,
      col='green3')

lines(RP.cluster.ul,hx,
      type='l',
      lty=2,
      lwd=1,
      col='green3')

lines(RP.mps,hx, # effect mps
      type='l',
      lwd=3,
      col='turquoise')

lines(RP.mps.ll,hx,
      type='l',
      lty=2,
      lwd=1,
      col='turquoise')

lines(RP.mps.ul,hx,
      type='l',
      lty=2,
      lwd=1,
      col='turquoise')

lines(RP.combined,hx, # combined effect
      type='l',
      lwd=3,
      col='deepskyblue3')

lines(RP.combined.ll,hx,
      type='l',
      lty=2,
      lwd=1,
      col='deepskyblue3')

lines(RP.combined.ul,hx,
      type='l',
      lty=2,
      lwd=1,
      col='deepskyblue3')


### end experiment 12 ###
########################


