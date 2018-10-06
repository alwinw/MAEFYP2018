#============================#
# Plotting Functions
# Alwin Wang
#----------------------------#

#--- Interpolation onto Regular Grid                              ----
#--- * Interpolate Single Variable                                ----
InterpGridVar <- function(redu, var, nx, ny, dp, linear, extrap) {
  # Rename column of interest
  colnames(redu)[colnames(redu) == var] = "z"
  # Interpolate onto regularly spaced grid
  grid <- interp(x = round(redu$x, dp),
                 y = round(redu$y, dp),
                 z = round(redu$z, dp),
                 nx = nx,
                 ny = ny,
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
dump_dump <- dump$dump; xvec = c(-0.5, 0.8); yvec = c(-0.2, 0.2); nx = 1500
dp = 12; linear = TRUE; extrap = FALSE
vars = c("u", "v", "p", "o")

InterpGrid <- function(dump_dump, vars, xvec, yvec, nx,
                       dp = 10, linear = TRUE, extrap = FALSE) {
  
  redu <- filter(dump_dump, x >= xvec[1], x <= xvec[2], y >= yvec[1], y <= yvec[2]) %>% 
    select(x, y, u, v, p, o) %>% 
    unique(.)
  ny = round(nx * (yvec[2] - yvec[1])/(xvec[2] - xvec[1]))
  
  grid_list <- lapply(vars, function(var) InterpGridVar(redu, var, nx, ny, dp, linear, extrap))
  grid   <- lapply(grid_list, function(df) df[3])
  grid   <- bind_cols(grid)
  grid$x <- grid_list[[1]]$x
  grid$y <- grid_list[[1]]$y
  
  ggplot()  +
    geom_raster (aes(x, y, fill = v),   grid,      interpolate = TRUE) +
    geom_polygon(aes(x, y),             dump$wall, fill = "white", colour = NA) +
    geom_path   (aes(x, y, colour = v), dump$wall, size = 0.2) +
    scale_fill_gradientn  (colours=rev(spectralpalette(10)), na.value = "white") +
    scale_colour_gradientn(colours=rev(spectralpalette(10)), na.value = "white",
                           limits = c(min(redu$v), max(redu$v)), guide = "none") +
    coord_fixed(xlim=c(-0.5, 0.8), ylim=c(-0.2, 0.2), expand = FALSE)
  
  ggplot()  +
    # geom_raster(aes(x, y, fill = u), grid, interpolate = TRUE) +
    # stat_contour(aes(x, y, z = u, fill = ..level..), filter(grid, !is.na(u)), 
                 # bins = 20, geom="polygon") +
    # geom_contour(aes(x, y, z = u), grid, bins = 20, colour = "white") +
    # geom_polygon(aes(x, y), dump$wall, fill = "white", colour = NA) +
    # geom_path(aes(x, y, colour = u), dump$wall, size = 0.2) +
    scale_fill_gradientn(
      colours=rev(spectralpalette(10)), na.value = "white") +
    scale_colour_gradientn(
      colours=rev(spectralpalette(10)), na.value = "white",
      limits = c(min(redu$u), max(redu$u)), guide = "none") +
    coord_fixed(xlim=c(-0.5, 0.8), ylim=c(-0.2, 0.2), expand = FALSE)
 
   
}



#--- Airfoil Plots                                                ----
