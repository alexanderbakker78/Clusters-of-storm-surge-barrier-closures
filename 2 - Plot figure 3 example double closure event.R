##########################################################################
#
# Plot example of double closure event with stochastic storm event (figure 3)
# Copyright (C) 2025  Alexander Bakker (TU Delft/Rijkswaterstaat)
#
# This script is developed as part of the study
#  "Storm Surge Clusters, Multi-Peaked Storms and Their Effect on
#   the Performance of the Maeslant Storm Surge Barrier (the  Netherlands)." 
#   Bakker, A.M.R.; Rovers, D.L.T.; Mooyaart, L.F. (2025)
#   J. Mar. Sci. Eng. 2025,13, x. https://doi.org/10.3390/xxxxx
#
# This script creates figure 3
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

rm(list=ls())

### load packages
require(plotrix)
require(greekLetters)

### input variables
slr       <- 0
hin_start <- 100               # [cmMSL]
A_basin   <- 152 * 1000 * 1000 # [m2]
hxx       <- 284
itide     <- 28
phi       <- -4
T         <- 71
Q         <- 2700

# tidal peak selection
file="Data/probability timing tidal peaks.csv"
tidalpeaks   <- read.table(file, header=TRUE, sep=";", dec=".")
tidalpeaks$X <- tidalpeaks$X+1 # index loopt van 0 t/m 697 (stom Python dingetje)
tidalpeaks$date_time <- as.POSIXct(tidalpeaks$date_time, format = "%Y-%m-%d %H:%M:%S")
tidalpeaks$tide_peak <- tidalpeaks$tide_peak + slr

# astronomic tide
file="Data/tides 2017.csv"
tides      <- read.table(file, header=TRUE, sep=",", dec=".")
tides$tide <- tides$tide + slr
tides$dt   <- as.POSIXct(tides$dt, format = "%Y-%m-%d %H:%M:%S")

tptide  <- tidalpeaks$date_time[itide]
tpsurge <- tptide + (phi*60*60)  
tphphi  <- tptide + (0.5*phi*60*60)

i=1

ns <- T*6 +1    # number of 10min timesteps
ts <- 1:ns      # time steps
  
  ## surge series ##
    surge <- hxx[i] * cos(pi * seq(-T[i]/2,T[i]/2,1/6) / T[i])^2
    
  ## tidal motion during surge ##
    tide1 <- which(tides$dt == tpsurge[i]) - (ns-1)/2
    tiden <- which(tides$dt == tpsurge[i]) + (ns-1)/2
    tide  <- tides$tide[tide1:tiden]
    time  <- tides$dt[tide1:tiden]
    
  ## water level
    swl   <- surge + tide
    hin   <- swl
    
  ## first closure
    if (max(swl) < 300) {  # no closure
      hx1 <- max(swl)
    } else {               # first closure
      exceedance1 <- which(swl >= 300)[1]
      cl1         <- tail(which(swl[1:exceedance1] < 200),1) + 1 # timing
      hin[cl1:ns] <- hin_start +
        (Q[i]/A_basin) * 10 * 60 * 100 * (ts[cl1:ns]-cl1)
      op1         <- (which(hin[cl1:ns]-swl[cl1:ns] >= 0) + cl1)[1]
      hx1      <- max(swl[cl1:op1-1])
      hin[op1:ns] <- swl[op1:ns]
      
      if (max(swl[op1:ns]) < 300) { # no second closure
        hx2 = max(swl[op1:ns])  
      } else {                      # second closure
        exceedance2 <- which(swl[op1:ns] >= 300)[1] + op1
        cl2         <- tail(which(swl[op1:exceedance2] < 200),1) + op1 + 1 # timing
        hin[cl2:ns] <- hin[cl2-1] +
          (Q[i]/A_basin) * 10 * 60 * 100 * (ts[cl2:ns]-cl2)
        op2         <- (which(hin[cl2:ns]-swl[cl2:ns] >= 0) + cl2)[1]
        hx2      <- max(swl[cl2:op2-1])
        hin[op2:ns] <- swl[op2:ns]
      }
    }
    
    midden <- length(surge+1)/2
    phase  <- midden + phi
    
    phase_diff <- time[c(phase,midden)]
    
    idtpsurge <- which(tides$dt == tpsurge)
    idtptide  <- which(tides$dt == tptide)
    
    plot(time,swl,type='l',xlab='date', ylab='Water level [cm MSL]', cex.lab=1.2)
    rect(tpsurge,-100,tptide,450,col = rgb(0.6,0.6,0.6,1/4), border=NA)
    lines(rep(tpsurge,2),c(-100,500), lty=2)
    lines(rep(tptide,2),c(-100,500), lty=2)
    arrows(x0=tpsurge, y0=150, x1=tptide, y1=150,
           length=0.1,
           angle=25,
           code=1)
    text(x=tphphi, y=155, pos=1, cex=1.0, paste(greeks("phi"),'= -4h'))

    lines(time,hin,type='l',col='red', lwd=2)
    lines(time,surge,type='l',col='green',lwd=2)
    lines(time,tide,type='l',col='blue',lwd=2)
    lines(time,swl,type='l',lwd=2)
    lines(time,rep(200,ns), lty=2 )
    lines(time,rep(300,ns), lty=2 )
    #text(x=time[10], y=215, pos=4, cex=1.1, "closure level (MSL +200cm)")
    #text(x=time[10], y=315, pos=4, cex=1.1, "closure decision level (MSL +300cm)")
    boxed.labels(time[50], 215, "closure level (MSL +200cm)", cex = 1.1, 
                 border = NA, bg ="white", xpad = 1.0, ypad = 1.4)
    boxed.labels(time[65], 315, "closure decision level (MSL +300cm)", cex = 1.1, 
                 border = NA, bg ="white", xpad = 1.0, ypad = 1.4)
    
    
        
    legend("topright", inset = 0.02,
           legend = c("Sea Water level", "Surge", "Tide", "Inner water level"),
           bg = 'white', # Legend background color
           box.col = "white",
           cex = 1.2,
           lty = c(1,1,1,1),
           col = c('black','green','blue','red'),
           lwd = 2)
    
    