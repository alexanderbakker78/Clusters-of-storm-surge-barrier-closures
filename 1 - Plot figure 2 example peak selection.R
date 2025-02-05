##########################################################################
#
# Plot example of peak selection (figure 2)
# Copyright (C) 2025  Alexander Bakker (TU Delft/Rijkswaterstaat)
#
# This script is developed as part of the study
#  "Storm Surge Clusters, Multi-Peaked Storms and Their Effect on
#   the Performance of the Maeslant Storm Surge Barrier (the  Netherlands)." 
#   Bakker, A.M.R.; Rovers, D.L.T.; Mooyaart, L.F. (2025)
#   J. Mar. Sci. Eng. 2025,13, x. https://doi.org/10.3390/xxxxx
#
# This script creates figure 2
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




##########################################################################
#
# 
# 
# Copyright (C) 2025  Alexander Bakker (TU Delft/Rijkswaterstaat)
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

# get data
file="Data/Drielingstorm HVH.csv"
data   <- read.table(file, header=TRUE, sep=";", dec=",")[,3:4]
swl    <- data[,2]
time   <- as.POSIXct(data[,1], format = "%d-%m-%Y %H:%M")
ns     <- length(time)

# plot data
plot(time,swl,type='l', col='blue', xlab='date', ylab='Water level [cm MSL]', lwd=2, cex.lab=1.2)

rect(time[140],-200,time[270],450,col = rgb(0.6,0.6,0.6,1/4), border=NA)
rect(time[360],-200,time[560],450,col = rgb(0.6,0.6,0.6,1/4), border=NA)
rect(time[725],-200,time[860],450,col = rgb(0.6,0.6,0.6,1/4), border=NA)

text(x=time[205], y=-90, pos=1, cex=1.3, "Dudley")
text(x=time[460], y=-90, pos=1, cex=1.3, "Eunice")
text(x=time[798], y=-90, pos=1, cex=1.3, "Franklin")

lines(time,rep(212.2,ns), lty=2 )
lines(time,rep(162.2,ns), lty=2 )
text(x=time[10], y=222.2, pos=4, cex=1.0, "threshold")
text(x=time[10], y=172.2, pos=4, cex=1.0, "end threshold")

arrows(x0=time[380], y0=240, x1=time[450], y1=223, lwd=2, col='red',
       length=0.15,
       angle=25,
       code=2)
arrows(x0=time[670], y0=255, x1=time[740], y1=262, lwd=2, col='red',
       length=0.15,
       angle=25,
       code=2)



