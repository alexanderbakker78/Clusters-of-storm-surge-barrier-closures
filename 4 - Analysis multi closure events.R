##########################################################################
#
# Analysis of multi closure events
# Copyright (C) 2025  Alexander Bakker (TU Delft/Rijkswaterstaat)
#
# This script is developed as part of the study
#  "Storm Surge Clusters, Multi-Peaked Storms and Their Effect on
#   the Performance of the Maeslant Storm Surge Barrier (the  Netherlands)." 
#   Bakker, A.M.R.; Rovers, D.L.T.; Mooyaart, L.F. (2025)
#   J. Mar. Sci. Eng. 2025,13, x. https://doi.org/10.3390/xxxxx
#
# This script performs a Monte Carlo simulation with a stochastic storm tide 
# event to explore the statistics of multi closure events as described in
# subsections 3.2 and 4.2
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

# load packages
#require(evmix)
#require(reshape2)
#require(extRemes)

# run supporting scripts with R-functions
source("4a - MPS_generate_random_input_MC_simulation.R")
source("4b - MPS_preselection_and_preperation_samples.R")
source("4c - MPS_simulate_storm.R")

# input/output files
ifile <- "Data/random_samples_applied.csv"
ofile <- "Results/random_samples.csv"

# input variables
nsamples  <- 524000
A_basin   <- 152 * 1000 * 1000  # [m2]
cdl       <- 300                # closure decision level
cl        <- 200                # closure level
slrs      <- c(0,25,50)         # sea level rise
hinstarts <- c(100,125,150,200) # initial average inner water level after closing

new_sample = F                  # T = sample new set, F = presampled set

###################################################
#   Monte Carlo Simulation                        #
###################################################
# random input samples
if(new_sample == T) {
  rsamples  <- generate_input(nsamples=nsamples)
  write.csv(rsamples, ofile, row.names=F)
} else {
  rsamples  <- read.csv(ifile, header=T, sep=',')
}

# output
freq <- data.frame()
vals <- seq(300,700,1)
cp.single <- data.frame(hx=vals) # conditional probabilities
cp.double <- data.frame(hx=vals) # conditional probabilities
cp.first  <- data.frame(hx=vals) # conditional probabilities
cp.second <- data.frame(hx=vals) # conditional probabilities
cp.third  <- data.frame(hx=vals) # conditional probabilities

# experiments
print('experiment')
experiment <- 0

for (slr in slrs) {
  samples   <- preselection(df=rsamples, slr=slr, cdl=cdl)
  
  for (hin_start in hinstarts) {
    experiment = experiment + 1
    expid      = paste('exp',as.character(experiment),sep='')
    print(experiment)
    
    # actual simulation
    df <- simulate_storm(df=samples, slr=slr, hin_start=hin_start, A_basin=A_basin, cl=cl, cdl=cdl)

    # sort
    single_closures <- df$hx1[which(is.na(df$hx2) & df$hx1 >= cdl)]
    first_closures  <- df$hx1[which(!is.na(df$hx2))]
    second_closures <- df$hx2[which(!is.na(df$hx2))]
    double_closures <- apply(data.frame(first_closures,second_closures),1,max)
    triple_closures <- df$hx3[which(!is.na(df$hx3))]
    
    # analysis
    freq <- rbind(freq, data.frame(exp=experiment,
                                   slr=slr,
                                   hin0=hin_start,
                                   single=length(single_closures),
                                   double=length(double_closures) -
                                     length(triple_closures),
                                   triple=length(triple_closures)))
    
    cp.single[expid] <- 1-ecdf(single_closures)(vals)
    cp.double[expid] <- 1-ecdf(double_closures)(vals)
    cp.first[expid]  <- 1-ecdf(first_closures)(vals)
    cp.second[expid] <- 1-ecdf(second_closures)(vals)
  }
}

freq[,4:6]             <- freq[,4:6]/100000
freq$total             <- freq$single + freq$double + freq$triple
freq$proportion_double <- freq$double/freq$total


## write output to csv
#write.csv(freq,"Results/frequency_closure_events.csv",row.names=F)
#write.csv(cp.single,"Results/cond_exc_prob_single.csv",row.names=F)
#write.csv(cp.double,"Results/cond_exc_prob_double.csv",row.names=F)
#write.csv(cp.first,"Results/cond_exc_prob_first.csv",row.names=F)
#write.csv(cp.second,"Results/cond_exc_prob_second.csv",row.names=F)

## figures
plot(cp.single$exp1, cp.single$hx, type='l', col='red2', lwd=3,
     cex.lab=1.4,
     cex.axis=1.4,
     log="x",
     xlim=c(1,1e-5),
     xlab='Conditonal exceedance probability',ylab='Sea Water Level [cmMSL]',
     xaxt = "n")
lines(cp.double$exp1, cp.single$hx, type='l', col='blue', lwd=5)
lines(cp.first$exp1,  cp.single$hx, type='l', col='green2', lwd=3)
lines(cp.second$exp1, cp.single$hx, type='l', col='orange', lwd=3)

lines(cp.single$exp4, cp.single$hx, type='l', col='red2', lwd=3,, lty=2)
lines(cp.double$exp4, cp.single$hx, type='l', col='blue', lwd=3, lty=2)
lines(cp.first$exp4,  cp.single$hx, type='l', col='green2', lwd=3, lty=2)
lines(cp.second$exp4, cp.single$hx, type='l', col='orange', lwd=3, lty=2)

lines(c(1,1e-5),c(360,360), lty=2 )
text(x=6.E-4, y=370, pos=4, cex=1.2, "critical interior water level [MSL 360cm]")

legend("topleft", inset = 0.01,
       legend = c("MSL +100cm                                    ",
                  "MSL +200cm                                    ",
                  "",
                  "Single closures","Double closures","First closures","Second closures"),
       title = "Initial water level behind barrier after closing",
      bg = 'white', # Legend background color
       box.col = "white",
       lty = c(1,2,1,1,1,1,1),
       col = c('black','black','white','red2','blue','green2','orange'),
       lwd = 2,
       cex=1.2)
axis(1, at = c(1,0.1,0.01,0.001,0.0001,0.00001),
     labels = c("1",
                expression("10"^{-1}),
                expression("10"^{-2}),
                expression("10"^{-3}),
                expression("10"^{-4}),
                expression("10"^{-5})),
     lwd = 0, lwd.ticks = 1.5, cex.axis=1.4)



