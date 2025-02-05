##########################################################################
#
# Sampling random input for Monte Carlo Simulation
# Copyright (C) 2025  Alexander Bakker (TU Delft/Rijkswaterstaat)
#
# This script is developed as part of the study
#  "Storm Surge Clusters, Multi-Peaked Storms and Their Effect on
#   the Performance of the Maeslant Storm Surge Barrier (the  Netherlands)." 
#   Bakker, A.M.R.; Rovers, D.L.T.; Mooyaart, L.F. (2025)
#   J. Mar. Sci. Eng. 2025,13, x. https://doi.org/10.3390/xxxxx
#
# This script generates the required number of random samples for the
# Monte Carlo simulation with a stochastic storm tide event as described in
# subsections 3.2
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

generate_input <- function(
    nsamples    = 524000,     # number of samples
    
    # input files
    file = "Data/probability timing tidal peaks.csv",
    
    # storm duration (lognormal distribution)
    meanT        = 54.3,        # mean storm duration
    varT         = 18.8^2,      # variance storm duration

    # surge height (generalize pareto distribution GPD)
    xi_hmx       =  -0.0062,
    sigma_hmx    =  31.456,
    th_hmx       =  107.6,
    
    # phase difference phi
    phase_values = -6.0:6.0,
    phase_probs  = c(0.006849, 0.027397, 0.054795, 0.095890,
                     0.102740, 0.171233, 0.034247, 0.041096,
                     0.034247, 0.054795, 0.123288, 0.164384, 0.089041),
    
    # river discharge Q (lognormal distribution)
    muQ          = 7.7,
    sigmaQ       = 0.5)

{
  # estimate parameters logn distribution storm duration
  muT     <- log(meanT / sqrt( (varT/meanT^2) + 1) ) #location parameter
  sigmaT  <- sqrt( log((varT/meanT^2) + 1)  )        # scale parameter
  
  # tidal peak selection
  tidalpeaks   <- read.table(file, header=TRUE, sep=";", dec=".")
  tidalpeaks$X <- tidalpeaks$X+1 # index loopt van 0 t/m 697 (stom Python dingetje)

  # sample surge height hx and tide number tide_i ###
  hxx    <- rgpd(n = nsamples, u = th_hmx, sigmau = sigma_hmx, xi = xi_hmx, phiu = 1)
  itide  <- sample(x = tidalpeaks$X, nsamples, replace = T, prob = tidalpeaks$probability)
  phi <- sample(x = as.numeric(phase_values), nsamples, replace = T, prob = phase_probs)
  
  T       <- round(rlnorm(n = nsamples, meanlog = muT, sdlog = sigmaT))
  Q       <- rlnorm(n = nsamples, meanlog = muQ, sdlog = sigmaQ)
  
  return(data.frame(hxx=hxx, itide=itide, phi=phi, T=T, Q=Q))
  
}

