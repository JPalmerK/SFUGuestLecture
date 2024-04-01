rm(list=ls())
library(argoFloats) # Load and manipulate the argos data
library(PAMmisc)
library(ocedata) # Required to plot argos data
library(ggplot2) # tidy package for plotting
library(scales) # also scales
library(viridisLite) # for colorblined friendly color scales

# Load the custom functions
source("~/GitHub/SFUOceanography/TL function.R")

# Load the argos data
load("~/GitHub/indexAll.RData")


plot(argosClean, which="map", bathymetry=TRUE,
     xlim=c(-129,-122), ylim=c(45,52),
     col=rgb(0,0,0,.5), cex=0.3, pch=16)

argosID = argosClean[['ID']]
index2 = subset(argosClean, ID="4901570") # select one track ID

pressureJan = index2[['argos']][[50]][['pressure']]
temperatureJan = index2[['argos']][[50]][['temperature']]
salJan = index2[['argos']][[50]][['salinity']]
par(mfrow=c(1,2))

pressureJul = index2[['argos']][[1]][['pressure']]
temperatureJul = index2[['argos']][[1]][['temperature']]
salJul = index2[['argos']][[1]][['salinity']]

# plot July profile
plot(temperatureJul, -pressureJul, ylim=rev(c(0,-1200)),
     xlim=c(3,13),
     xlab='Temperature (C)', ylab='Depth (dbar)',
     pch=16, col=rgb(0,0,0,.4),
     type='o')
legend('bottomright', legend='July',bty='n')


# plot Jan profile
plot(temperatureJan, -pressureJan, ylim=rev(c(0,-1200)),
     xlim=c(3,13),
     xlab='Temperature (C)', ylab='Depth (dbar)',
     pch=16, col=rgb(0,0,0,.4),
     type='o')
legend('bottomright', legend='January',bty='n')


# Equation for soundspeed based on temperature, salinity and depth
soundSpeedMacenzie<-function(D,S,Temp){
  c = 1448.96 + 4.591*Temp - 5.304 * 10^-2*Temp^2 + 2.374 * 10^-4*Temp^3 +
      1.340*(S-35) + 1.630 * 10^-2*D + 1.675 * 10^-7*D^2 - 1.025 *
    10^-2*Temp*(S - 35)}

# Functions to create a smooth temp/salintiy gradient.
salJulfx = approxfun(pressureJul, salJul)
tempJulfx = approxfun(pressureJul, temperatureJul)


# Set up the soundspeed profile
depth <- seq(0, 1200, by =3)
temp <- tempJulfx(depth)
sal <- salJulfx(depth)

# Test by plotting
plot(temp, -depth)
plot(sal, -depth)


# Use the salinity and temperature to create a soundspeed profile
ssp <-soundSpeedMacenzie(D = depth, S = sal, Temp = temp)

# check it again
plot(ssp, -depth)

# the argos float only goes down so deep (hitting the bottom with an expensive
# bit of kit is ill advised). So lets populate the rest of the SSP with the value
# of the last observed soundspeed. This isn't perfect but it's a resonable 
# approximation

# Clear out SSP values with NA
depth= depth[!is.na(ssp)]; ssp=ssp[!is.na(ssp)]

rt<-raytrace_TL(x0 = 0, z0 = 5, theta0 = 2, tt = 60, 
             zz = depth, cc = ssp, TRUE)


# Create multiple rays- use the sequence command to creat an array from -pi to 
# pi
theta <-seq(-90, 90, length.out=20)


# Generate data using the raytrace function
data <- raytrace_TL(x0 = 0, z0 = 200, theta0 = theta, tt = 60, 
                    zz = depth, cc = ssp, plot = FALSE)


min_length <- min(lengths(data$x), lengths(data$z), lengths(data$TL))

# Truncate x, z, and TL vectors to the minimum length
x <- unlist(lapply(data$x, head, n = min_length))
z <- unlist(lapply(data$z, head, n = min_length))
TL <- unlist(lapply(data$TL, head, n = min_length))


dfout= data.frame(x=x, z =z, TL=TL)
dfout= subset(dfout, x<2000)


# Plot using ggplot and geom_point with color gradient
ggplot(dfout, aes(x = x, y = -z, color = TL)) +
  geom_point() +
  scale_color_viridis_c()+  # Adjust color scale as needed
  labs(x = "Range (m)", y = "Depth (m)", color = "Average TL") +
  ggtitle("Average Transmission Loss (TL) by Depth and Range")


# The above code shows the rays for several anglesand the reflections and 
# refractions. To estimate the overall TL we need to use the acoustic recriprosity
# principal which states that source and receiver locations can be switched 
# without impacting the received level of the sound. The following function
# will take in our grid space that we want and estimate the overall transmission
# loss between the source at the initial sensor location and a hypothetical receiver
# at the each grid cell


# Example usage with user-defined parameters
soureX <- 0  # Source x-coordinate
sourceDepth <- 20  # Source depth
theta_range <- seq(-90, 90, by = 1)  # Range of angles for multiple rays
tt <- 100  # Total travel time

# Set up the intended output grid for the transmission loss
range_extent <- 20000  # Range extent in meters (e.g., 20 km)
depth_extent <- max(depth)  # Depth extent in meters
range_resolution <- 100  # Range resolution in meters (e.g., 100m)
depth_resolution <- 10  # Depth resolution in meters (e.g., 10m)

# Run the transmission loss model, this will take a little while
result_grid <- compute_avg_TL_grid(soureX, 
                                   sourceDepth, 
                                   theta_range, tt, depth, ssp,
                                   range_extent, depth_extent,
                                   range_resolution, depth_resolution,
                                   progress = TRUE,
                                   plot = TRUE)


# Plot the result
p<- plotTLgrid(result_grid)
p