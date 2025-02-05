##########################################################################
#
# Analysis of multi closure events
# Copyright (C) 2025  Alexander Bakker (TU Delft/Rijkswaterstaat)
#
# This function is developed as part of the study
#  "Storm Surge Clusters, Multi-Peaked Storms and Their Effect on
#   the Performance of the Maeslant Storm Surge Barrier (the  Netherlands)." 
#   Bakker, A.M.R.; Rovers, D.L.T.; Mooyaart, L.F. (2025)
#   J. Mar. Sci. Eng. 2025,13, x. https://doi.org/10.3390/xxxxx
#
# This functions performs a deterministic simulation with the 'stochastic'
# storm tide event as described in section 3.2
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

simulate_storm <- function(
    
    #random samples
    df, 

    # input files
    file = "Data/tides 2017.csv",

    # input variables  
    slr       = 0,                 # sea level rise [cm]
    hin_start = 100,               # initial inner water level after closing [cmMSL]
    A_basin   = 300 * 1000 * 1000, # surface area inner basin [m2]
    cdl       = 300,               # closure decision level
    cl        = 200)               # closure level

{

  # extract samples
  nsamples <- nrow(df)
  hxx      <- df$hxx
  itide    <- df$itide
  phi      <- df$phi
  tptide   <- df$tptide
  tpsurge  <- df$tpsurge
  T        <- df$T
  Q        <- df$Q
  
  # astronomic tide
  file="Data/tides 2017.csv"
  tides      <- read.table(file, header=TRUE, sep=",", dec=".")
  tides$tide <- tides$tide + slr
  tides$dt   <- as.POSIXct(tides$dt, format = "%Y-%m-%d %H:%M:%S")
  
  # outputs
  hx1 <- rep(NA,nsamples)
  hx2 <- rep(NA,nsamples)
  hx3 <- rep(NA,nsamples)
  
  # simulate samples
  for (i in 1:nsamples) {
    ns <- T[i]*6 +1 # number of 10min time steps
    ts <- 1:ns      # time steps
    
    ## surge series ##
    surge <- hxx[i] * cos(pi * seq(-T[i]/2,T[i]/2,1/6) / T[i])^2
    
    ## tidal motion during surge ##
    tide1 <- which(tides$dt == tpsurge[i]) - (ns-1)/2
    tiden <- which(tides$dt == tpsurge[i]) + (ns-1)/2
    tide  <- tides$tide[tide1:tiden]
    
    ## water level
    swl   <- surge + tide
    hin   <- swl
    
    ## first closure
    if (max(swl) < cdl) {                                         # no closure
      hx1[i] <- max(swl)                                            # max wl in case of no closure
    } else {               # first closure                        # first/single closure
      exceedance1 <- which(swl >= cdl)[1]                           # timing 1st exceedance cdl
      cl1         <- tail(which(swl[1:exceedance1] < cl),1) + 1     # timing 1st closure
      hin[cl1:ns] <- hin_start +
        (Q[i]/A_basin) * 10 * 60 * 100 * (ts[cl1:ns]-cl1)           # inner water level 1st cl
      op1         <- which(hin[cl1:ns] >= swl[cl1:ns])[1] + cl1 - 1 # timing of 1st opening

      if(is.na(op1) | op1 > ns) {                                   # no opening within surge duration
        hx1[i] <- max(swl)                                          # max water level 1st closure
      } else {
        hx1[i]      <- max(swl[cl1:op1-1])                          # max water level 1st closure
        hin[op1:ns] <- swl[op1:ns]                                  

        ## second closure
        if (max(swl[op1:ns]) >= cdl) {                             # second closure
          exceedance2 <- which(swl[op1:ns] >= cdl)[1] + op1 - 1         # timing 2nd exceedance cdl
          cl2         <- tail(which(swl[op1:exceedance2] < cl),1) + op1 # timing 2nd closure
          if(length(cl2) == 0) cl2 <- op1 + 1                           # if swl at opening exceeds cl
          
          hin[cl2:ns] <- swl[cl2 - 1] +
            (Q[i]/A_basin) * 10 * 60 * 100 * (ts[cl2:ns]-cl2)           # inner water level 2nd cl
          op2         <- which(hin[cl2:ns] >= swl[cl2:ns])[1] + cl2 -1  # timing of 2nd opening

          if(is.na(op2) | op2 > ns) {                                   # no opening within surge duration
            hx2[i] <- max(swl)                                          # max water level 2nd closure
          } else {
            hx2[i]      <- max(swl[cl2:op2-1])                          # max water level 2nd closure
            hin[op2:ns] <- swl[op2:ns]                                  
            
            ## thrid closure (not further elaborated)
            if (max(swl[op2:ns]) >= cdl) hx3[i] = max(swl[op2:ns])      # max water level 3rd closure
          } # endif 2nd opening
        } # endif 2nd closure
        
      } # endif 1st opening
    } # endif 1st/only closure
  } # for (i in 1:nsamples)
  
  return(cbind(slr,hin_start,df,hx1,hx2,hx3))
  
}



