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
# This functions performs a preselection of the random samples to reduce 
# the computational burden in the Monte Carlo simulation
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

preselection <- function(
    df,
#    nsamples    = 524000,     # number of samples
    
    # input files
    file = "Data/probability timing tidal peaks.csv",
    
    # input variables
    cdl = 300,                  # closure decision level
    slr = 0)                    # sea level rise
    
{
  # tidal peak selection
  tidalpeaks   <- read.table(file, header=TRUE, sep=";", dec=".")
  tidalpeaks$X <- tidalpeaks$X+1 # index loopt van 0 t/m 697 (stom Python dingetje)
  tidalpeaks$date_time <- as.POSIXct(tidalpeaks$date_time, format = "%Y-%m-%d %H:%M:%S")
  tidalpeaks$tide_peak <- tidalpeaks$tide_peak + slr
  
  # preselect samples that might require a closure
  index  <- which(df$hxx + tidalpeaks$tide_peak[df$itide] >= cdl)
  
  ns      <- length(index)         # number of samples that require further analysis
  hxx     <- df$hxx[index]
  itide   <- df$itide[index]
  phi     <- df$phi[index]
  T       <- df$T[index]
  Q       <- df$Q[index]
  
  tptide  <- tidalpeaks$date_time[itide] # timing tidal peak
  tpsurge <- tptide + (phi*60*60)        # timing surge peak
  
  return(data.frame(hxx=hxx, itide=itide, phi=phi, tptide=tptide, tpsurge=tpsurge, T=T, Q=Q))
}

