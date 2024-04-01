library(viridisLite)
library(ggplot2)




raytrace_TL <- function(x0, z0, theta0, tt, zz, cc, plot = TRUE, progress = FALSE) {
  # Raytracing function based on Taiki Sakai's R implementation
  
  # Add a zero depth and corresponding sound speed if not already present
  if (zz[1] != 0) {
    zz <- c(0, zz)
    cc <- c(cc[1], cc)
  }
  
  # Initialize variables and constants
  MAX_BOUNCES <- 20
  n_theta <- length(theta0)
  ray_paths_x <- vector('list', length = n_theta)
  ray_paths_z <- ray_paths_x
  ray_paths_d <- ray_paths_x
  ray_paths_t <- ray_paths_x
  ray_paths_tt <- ray_paths_x
  ray_paths_theta <- ray_paths_x
  transmission_loss <- ray_paths_x
  
  # Iterate over each angle of incidence
  for (m in seq_along(theta0)) {
    n_svp <- length(zz)
    zz_end <- zz[n_svp]
    
    dz <- diff(zz)
    dc <- diff(cc)
    
    # Check if the ray is coming from above or below
    if (theta0[m] < 0 && z0 > 0) {
      z <- rev(zz)
      zzz <- c(zz[1], zz[1] + cumsum(rev(dz)))
      ccc <- c(cc[n_svp], cc[n_svp] - cumsum(rev(dc)))
      mult <- 1
    } else {
      z <- zz
      zzz <- c(zz[1], zz[1] + cumsum(dz))
      ccc <- c(cc[1], cc[1] + cumsum(dc))
      mult <- -1
    }
    dz0 <- dz
    
    # Extend ray path until it reaches MAX_BOUNCES times the maximum depth
    while (abs(zzz[length(zzz)]) < MAX_BOUNCES * zz_end) {
      if (mult == -1) {
        zzz <- c(zzz, zzz[length(zzz)] + cumsum(rev(dz)))
        ccc <- c(ccc, ccc[length(ccc)] + mult * cumsum(rev(dc)))
      } else {
        zzz <- c(zzz, zzz[length(zzz)] + cumsum(dz))
        ccc <- c(ccc, ccc[length(ccc)] + mult * cumsum(dc))
      }
      z <- c(z, z[length(z)] + mult * cumsum(rev(dz0)))
      dz <- rev(dz)
      dz0 <- rev(dz0)
      mult <- mult * -1
    }
    
    dz <- diff(zzz)
    g <- diff(ccc) / dz
    
    idx <- which.min(abs(z[1:n_svp] - z0))
    q <- idx[1]
    qstart <- q
    
    theta <- rep(NA, length(dz))
    theta[q] <- abs(theta0[m]) * pi / 180
    x <- rep(NA, n_svp)
    dx <- rep(NA, n_svp)
    x[q] <- x0
    
    t <- 0
    d <- 0
    TL <- rep(NA, n_svp)
    
    if (progress) {
      pb <- txtProgressBar(min = 0, max = tt, style = 3)
    }
    
    while (t < tt) {
      if (q > length(g)) {
        warning('CCOM:outofBounds Not enough bounces specified for time requested. Increase MAXBOUNCES')
        break
      }
      
      if (g[q] == 0) {
        theta[q+1] <- theta[q]
        if (theta[q] == 0) {
          dx[q] <- abs(dz[q])
        } else {
          dx[q] <- abs(dz[q] / tan(theta[q]))
        }
        dd <- sqrt(dx[q]^2 + dz[q]^2)
        dt <- dd / ccc[q]
      } else {
        Rc <- -1/g[q] * ccc[q] / cos(theta[q])
        tmpCos <- cos(theta[q]) - dz[q] / Rc
        
        if (abs(tmpCos) > 1) {
          theta[q+1] <- -theta[q]
          tmp <- tryCatch({
            which.min(abs(z[(q+1):(q+2*length(cc))] - z[q]))
          }, error = function(e) {
            numeric(0)
          })
          if (length(tmp) == 0) {
            warning('CCOM:outofBounds Not enough bounces specified for time requested. Increase MAXBOUNCES')
            break
          }
          
          ccc <- c(ccc[1:q], ccc[(q + tmp[1]):length(ccc)])
          zshift <- zzz[q + tmp[1] + 1] - zzz[q]
          zzz <- c(zzz[1:q], zzz[(q + tmp[1]):length(zzz)] - zshift)
          z <- c(z[1:q], z[(q + tmp[1]):length(z)])
          dz <- diff(zzz)
          g <- diff(ccc) / dz
          dx[q] <- 0
          dt <- 0
          dd <- 0
        } else {
          theta[q+1] <- acos(tmpCos)
          if (theta[q] > 0) {
            dx[q] <- Rc * (sin(theta[q+1]) - sin(theta[q]))
            dt <- -2/g[q] * (atanh(tan(theta[q+1]/2)) - atanh(tan(theta[q]/2)))
          } else {
            dx[q] <- abs(Rc * -1*(sin(theta[q+1]) - sin(theta[q])))
            dt <- -2/g[q] * (atanh(tan(theta[q+1]/2)) - atanh(tan(theta[q]/2)))
          }
          dd <- dt * (ccc[q] + g[q]/2)
        }
      }
      
      t <- t + dt
      x[q+1] <- x[q] + dx[q]
      d <- d + dd
      q <- q + 1
      TL[q] <- 18 * log10(d)
      if (progress) {
        setTxtProgressBar(pb, value = t)
      }
    }
    
    x <- x[!is.na(x)]
    z <- z[qstart:(qstart + length(x) - 1)]
    TL <- TL[!is.na(TL)]
    
    corrector <- (tt - (t - dt)) / dt
    x[length(x)] <- x[length(x) - 1] + (x[length(x)] - x[length(x) - 1]) * corrector
    z[length(z)] <- z[length(z) - 1] + (z[length(z)] - z[length(z) - 1]) * corrector
    TL[length(TL)] <- TL[length(TL) - 1] + (TL[length(TL)] - TL[length(TL) - 1]) * corrector
    
    t <- (t - dt) + dt * corrector
    theta <- theta[!is.na(theta)]
    ray_paths_x[[m]] <- x
    ray_paths_z[[m]] <- z
    ray_paths_d[[m]] <- d
    ray_paths_t[[m]] <- t
    ray_paths_theta[[m]] <- theta * 180 / pi
    transmission_loss[[m]] <- TL
  }
  
  # Plotting the results if requested
  if (length(plot) == 1) {
    plot <- rep(plot, 2)
  }
  
  if (plot[1]) {
    plot(cc[1:n_svp], zz, ylab = 'Depth, m', xlab = 'Sound Speed, m/s', type = 'l')
    title('Sound Speed Profile')
  }
  
  if (plot[2]) {
    x_range <- range(sapply(ray_paths_x, range))
    z_range <- range(sapply(ray_paths_z, range))
    cols <- viridis_pal()(32)
    col_scale <- 31 * ((theta0 - min(theta0)) / max(theta0 - min(theta0))) + 1
    if (length(theta0) == 1) {
      col_scale <- 1
    }
    
    for (i in seq_along(ray_paths_x)) {
      if (i == 1) {
        plot(ray_paths_x[[i]], -ray_paths_z[[i]], type = 'l', xlim = x_range, ylim = rev(-z_range), col = cols[col_scale[i]],
             xlab = 'Range (m)', ylab = 'Depth (m)', main = 'Ray Paths')
      } else {
        lines(ray_paths_x[[i]], -ray_paths_z[[i]], col = cols[col_scale[i]])
      }
    }
    
    if (length(theta0) > 1) {
      leg <- rep(NA, 11)
      leg[c(1, 6, 11)] <- c(min(theta0), mean(range(theta0)), max(theta0))
      legend(x = 'topright',
             legend = c(NA, leg),
             fill = c(NA, viridis_pal()(11)),
             border = NA,
             y.intersp = .5,
             cex = 1,
             text.font = 1,
             title = 'Start Angle'
      )
    }
  }
  
  # Return the results as a list
  list(x = ray_paths_x, z = ray_paths_z, t = ray_paths_t, d = ray_paths_d, TL = transmission_loss)
}



compute_avg_TL_grid <- function(x0, z0, theta_range, tt, zz, cc, range_extent, 
                                depth_extent, range_resolution, 
                                depth_resolution, progress = FALSE, plot = TRUE) {
  
  
  # Step 1: Generate ray trace paths for all theta values
  ray_paths <- list()
  for (theta in theta_range) {
    result <- raytrace_TL(x0, z0, theta, tt, zz, cc, plot = FALSE, progress = progress)
    ray_paths[[as.character(theta)]] <- result
  }
  
  # Step 2: Initialize transmission loss grid
  range_grid <- seq(0, range_extent, by = range_resolution)
  depth_grid <- seq(0, depth_extent, by = depth_resolution)
  avg_TL_grid <- matrix(NA, nrow = length(range_grid), ncol = length(depth_grid))
  
  # To count how many rays contribute to each grid cell
  ray_count_grid <- matrix(1, nrow = length(range_grid), ncol = length(depth_grid))    
  
  # Step 3: Populate transmission loss grid based on ray intersections
  for (theta in theta_range) {
    ray_path <- ray_paths[[as.character(theta)]]
    x_values <- ray_path$x[[1]]  # Extract as vector
    z_values <- ray_path$z[[1]]  # Extract as vector
    TL_values <- ray_path$TL[[1]]
    
    # Step through the x values of the range grid and get Y values
    for (ii in 1:length(range_grid)-1) {
      
      # all the x values in range of the first grid cell
      idxInrange = which(x_values>= range_grid[ii] & x_values<= range_grid[ii+1])
      
      # get the y and TL values associated with the range values
      yvalues = z_values[idxInrange]
      tlvalues <- TL_values[idxInrange]
      
      # get the index of the y (depth) value son the TL grid
      yGrididx = which(depth_grid>= min(yvalues) & depth_grid <= max(yvalues))
      
      
      # now step through the y grid indicies and populate the TL values
      for(kk in 1:(length(yGrididx)-1)){
        
        ygridMin = depth_grid[yGrididx[kk]]
        ygridMax = depth_grid[yGrididx[kk+1]]
        
        # get the TL values and convert to linear
        tlkeep = tlvalues[yvalues<= ygridMax & yvalues>= ygridMin]
        tlkeep = mean(10^(tlkeep/20), na.rm = TRUE)
        
        #Get the tL values that fell within the delta y
        avg_TL_grid[ii, yGrididx[kk]]=  tlkeep
        ray_count_grid[ii, yGrididx[kk]]= ray_count_grid[ii, yGrididx[kk]]+1
      }
      
      
    }
  }
  
  # Do the averages for overlapping TL values and convert back to dB
  dbTl= 20*log10(avg_TL_grid/ray_count_grid)
  
  # Return the transmission loss grid
  return(list(range = range_grid, depth = depth_grid, avg_TL = dbTl))
}


# Plot the transmission loss grid
plotTLgrid<-function(result_grid){
  
  # define jet colormap
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                                   "cyan", "#7FFF7F", "yellow", "#FF7F00",
                                   "red", "#7F0000"))
  
  # Convert avg_TL_grid to long format manually
  long_data <- expand.grid(Range = result_grid$range, 
                           Depth = -result_grid$depth)
  long_data$TransmissionLoss <- as.vector(result_grid$avg_TL)
  
  # Plot with ggplot using the long format data
  p<-ggplot(long_data, aes(x = Range, y = Depth, fill = TransmissionLoss)) +
    geom_tile() +
    scale_fill_gradientn(colors = jet.colors(7),
                         limits=c(20,80),
                         breaks= floor(seq(20,80, length.out = 4)))+
    labs(x = "Range (m)", y = "Depth (m)", 
         title = "Average Transmission Loss Grid",
         fill="Transmission Loss (dB)")
  
  return(p)

  
  breaks=c(0,0.5,1)
}

# Calculat the the speed of sound using Mcenzy
# Equation for soundspeed based on temperature, salinity and depth
soundSpeedMackenzie<-function(D,S,Temp){
  c = 1448.96 + 4.591*Temp - 5.304 * 10^-2*Temp^2 + 2.374 * 10^-4*Temp^3 +
    1.340*(S-35) + 1.630 * 10^-2*D + 1.675 * 10^-7*D^2 - 1.025 *
    10^-2*Temp*(S - 35)}

