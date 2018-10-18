#============================#
# Plotting Functions
# Alwin Wang
#----------------------------#

# ggplot2 setup (consider moving to a separate script)
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 

#--- Interpolation onto Regular Grid                              ----
#--- * Interpolate Single Variable                                ----
InterpGridVar <- function(redu, var, nx, ny, dp, linear, extrap) {
  # Rename column of interest
  colnames(redu)[colnames(redu) == var] = "z"
  # Interpolate onto regularly spaced grid
  grid <- interp(x = round(redu$x, dp),
                 y = round(redu$y, dp),
                 z = round(redu$z, dp),
                 nx = nx, ny = ny,
                 duplicate = "strip",
                 linear = linear,
                 extrap = extrap)
  # Converge to dataframe
  grid <- interp2xyz(grid, data.frame = TRUE)
  # Remove/round any outliers
  grid <- grid %>% 
    mutate(z = ifelse(z < min(redu$z), min(redu$z), z)) %>% 
    mutate(z = ifelse(z > max(redu$z), max(redu$z), z))
  # Rename the column
  colnames(grid)[3] <- var
  return(grid)
}

#--- * Interpolate All Variables onto Grid                        ----
# dump_dump <- dump$dump; xvec = c(-0.5, 0.8); yvec = c(-0.2, 0.2); nx = 1500
# dp = 12; linear = TRUE; extrap = FALSE
# vars = c("u", "v", "p", "o")
InterpGrid <- function(dump_dump, vars, xvec, yvec, nx,
                       dp = 10, linear = TRUE, extrap = FALSE) {
  # Increase domain to allow for NA at edges/corners
  xveci = xvec + c(-0.1*(xvec[2]-xvec[1]), 0.1*(xvec[2]-xvec[1]))
  yveci = yvec + c(-0.1*(yvec[2]-yvec[1]), 0.1*(yvec[2]-yvec[1]))
  # Reduce to range in xvec and yvec
  redu <- filter(dump_dump, x >= xveci[1], x <= xveci[2], y >= yveci[1], y <= yveci[2]) %>% 
    select(x, y, u, v, p, o) %>% 
    unique(.)
  # Number of points in y direction for aspect ratio = 1.5
  ny = round(nx * (yveci[2] - yveci[1])/(xveci[2] - xveci[1]) * 1.5)
  # Interpolate and clean up
  grid_list <- lapply(vars, function(var) InterpGridVar(redu, var, nx, ny, dp, linear, extrap))
  grid   <- lapply(grid_list, function(df) df[3])
  grid   <- bind_cols(grid)
  grid$x <- grid_list[[1]]$x
  grid$y <- grid_list[[1]]$y
  grid <- grid[c("x", "y", colnames(grid)[1:(ncol(grid)-2)])]
  # Return result
  return(grid)
}

#--- Airfoil Plots                                                ----
# dump_dump <- dump$dump; scalearr = c(5, 2);
PlotAirfoil <- function(dump_dump, dump_wall, data_plot, scalearr = c(10, 4), saveplot) {
  # Ranges
  xvec = c(-0.5, 0.8)
  yvec = c(-0.2, 0.2)
  vars = c("u", "v", "p", "o")
  # Interpolate
  grid <- InterpGrid(dump_dump, vars, xvec, yvec, nx = 500, dp = 12)
  grid <- filter(grid, !is.na(u), !is.na(v))
  # Create end points
  grid <- grid %>% 
    mutate(vlen = EucDist(u, v))
  # Scale the output for the arrow plotting to be different to interpgrid
  scale  = (grid$x[2] - grid$x[1])/max(grid$vlen)*scalearr[1]
  scalex = sort(unique(grid$x))
  scaley = sort(unique(grid$y))
  scalex = scalex[seq(1, length(scalex), scalearr[1])]
  scaley = scaley[seq(1, length(scaley), scalearr[2])]
  # Magnitude of the velocity
  plot_vel <- ggplot() +
    geom_raster(
      aes(x, y, fill = sqrt(u^2+v^2)), grid, 
      interpolate = TRUE) +
    scale_fill_gradientn(
      colours=rev(spectralpalette(10)), na.value = "white",
      limits = c(0, 1.5),
      name = expression(group("|",bar(u),"|"))) +
    coord_fixed(xlim=xvec, ylim=yvec, expand = FALSE) +
    geom_polygon(aes(x, y), dump_wall, fill = "white", colour = NA) +
    geom_path   (aes(x, y), dump_wall, colour = "black", size = 0.2) +
    ggtitle(expression(paste("Velocity,  ", group("|",bar(u),"|"), " = ", sqrt(u^2 + v^2))))
  # Vector field over vorticity
  plot_vor <- ggplot() +
    geom_raster(
      aes(x, y, fill = o), grid, 
      interpolate = TRUE) +
    scale_fill_gradientn(
      colours=rev(spectralpalette(10)), na.value = "white", 
      limits = c(-650, 600),
      name = expression(omega)) +
    geom_segment(
      aes(x, y, xend = x + u*scale, yend = y + v*scale), 
      filter(grid, x%in%scalex,  y%in%scaley),
      size = 0.2) +
    geom_polygon(aes(x, y), dump_wall, fill = "white", colour = NA) +
    geom_path   (aes(x, y), dump_wall, colour = "black", size = 0.2) +
    coord_fixed(xlim=xvec, ylim=yvec, expand = FALSE) +
    ggtitle("Vorticity")
  # Sizes
  # format(object.size(plot_vel), units = "auto")
  # format(object.size(plot_vor), units = "auto")
  # Save Images
  subsavefold <- paste0(saveplot, "af/")
  if(!dir.exists(subsavefold)) dir.create(subsavefold)
  save_path <- paste0(subsavefold, data_plot$plotname)
  ggsave(paste0(save_path, "_af-vel.png"), plot_vel,
         scale = 1.5, width = 10, height = 4, units = "cm", dpi = 600)
  ggsave(paste0(save_path, "_af-vor.png"), plot_vor,
         scale = 1.5, width = 10, height = 4, units = "cm", dpi = 600)
  # Return plot
  return(list(afvel = plot_vel, afvor = plot_vor))
}

# dump_dump <- dump$dump; dump_wall <- dump$wall
PlotTE <- function(dump_dump, dump_wall, data_plot, scalearr = c(10, 4), saveplot) {
  # Ranges
  xvec = c( 0.595,  0.615)
  yvec = c(-0.046, -0.030)
  vars = c("u", "v", "o")
  # Interpolate
  grid <- InterpGrid(dump_dump, vars, xvec, yvec, nx = 500, dp = 12, linear = TRUE)
  grid <- filter(grid, !is.na(u), !is.na(v))
  # Create end points
  grid <- grid %>% 
    mutate(vlen = EucDist(u, v))
  # Scale the output for the arrow plotting to be different to interpgrid
  scale  = (grid$x[2] - grid$x[1])/max(grid$vlen)*scalearr[1]
  scalex = sort(unique(grid$x))
  scaley = sort(unique(grid$y))
  scalex = scalex[seq(1, length(scalex), scalearr[1])]
  scaley = scaley[seq(1, length(scaley), scalearr[2])]
  # Vector field over vorticity
  plot_vor <- ggplot() +
    geom_raster(
      aes(x, y, fill = o), grid, 
      interpolate = TRUE) +
    scale_fill_gradientn(
      colours=rev(spectralpalette(10)), na.value = "white", 
      limits = c(-250, 250),
      name = expression(omega)) +
    geom_segment(
      aes(x, y, xend = x + u*scale, yend = y + v*scale), 
      filter(grid, x%in%scalex,  y%in%scaley),
      size = 0.2) +
    geom_polygon(aes(x, y), dump_wall, fill = "white", colour = NA) +
    geom_path   (aes(x, y), dump_wall, colour = "black", size = 0.2) +
    # geom_point  (aes(x, y, fill = o), dump_dump, shape = 21) +
    coord_fixed(xlim=xvec, ylim=yvec, expand = FALSE) +
    ggtitle("Vorticity")
  
  # Sizes
  # format(object.size(plot_vel), units = "auto")
  # format(object.size(plot_vor), units = "auto")
  # Save Images
  subsavefold <- paste0(saveplot, "te/")
  if(!dir.exists(subsavefold)) dir.create(subsavefold)
  save_path <- paste0(subsavefold, data_plot$plotname)
  ggsave(paste0(save_path, "_te-vor.png"), plot_vor,
         scale = 1.5, width = 10, height = 8, units = "cm", dpi = 300)
  # Return plot
  return(list(tevor = plot_vor))
}











