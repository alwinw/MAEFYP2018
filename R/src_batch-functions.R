#============================#
# Batch Functions
# Alwin Wang
#----------------------------#

#--- Airfoil Calculation                                          ----
# Determine things like spline distance once per unique airfoil i.e. boundary profile and wall output
BatchLoadAirfoil <- function(data_airfoil, srcpath = "") {
  source(paste0(srcpath, "src_library-manager.R"))                # Call libraries and install missing ones
  source(paste0(srcpath, "src_helper-functions.R"))               # Smaller functions used
  #--- * Boundary Data                                              ----
  # Not sure this is actually used anywhere...
  bndry <- LoadBndry(data_airfoil$folder)
  #--- * Wall Mesh Data                                             ----
  wallmesh <- LoadWallGrad(data_airfoil$seshpath)
  long_wall <- AirfoilLongWall(wallmesh)
  long_wall <- AirfoilSpline(long_wall)
  long_wall <- AirfoilNorm  (long_wall)
  #--- > Airfoil Calc Output                                        ----
  list_airfoil <- list(
    airfoil   = data_airfoil,
    bndry     = bndry,
    long_wall = long_wall)
  rm(data_airfoil, bndry, long_wall, wallmesh)
  # Output
  return(list_airfoil)
}

#--- Session and Mesh Calculation                                 ----
BatchLoadMesh <- function(data_mesh, outp_airfoil, srcpath = "") {
  source(paste0(srcpath, "src_library-manager.R"))                # Call libraries and install missing ones
  source(paste0(srcpath, "src_helper-functions.R"))               # Smaller functions used
  list_airfoil <- outp_airfoil[[data_mesh$airfoil]]
  #--- * Airfoil Data                                               ----
  long <- list()
  #--- * Session Data                                               ----
  session   <- LoadSeshFileKeywords(data_mesh$seshpath)
  long$sesh <- LongSesh(session)
  #--- * Mesh Data                                                  ----
  long$mesh <- LoadMesh(data_mesh$seshpath)
  long$mesh <- LongMesh(long$mesh, long$sesh)
  #--- * Wall Data                                                  ----
  long$wall <- list_airfoil$long_wall
  long$wall <- LongWall(long$wall, long$mesh)
  #--- * Local Data                                                 ----
  long$mesh <- LocalMesh(long$mesh, long$wall)
  #--- > Sesh & Mesh Calc Output                                    ----
  list_mesh <- list(
    wall = long$wall,
    mesh = long$mesh)
  rm(data_mesh, long, session)
  # output
  return(list_mesh)
}

#--- Dump File Calculation                                        ----
BatchLoadDump <- function(data_dump, outp_mesh, plot = 0, srcpath = "") {
  source(paste0(srcpath, "src_library-manager.R"))                # Call libraries and install missing ones
  source(paste0(srcpath, "src_helper-functions.R"))               # Smaller functions used
  list_mesh <- outp_mesh[[data_dump$ID]]
  #--- * Dump Data                                                  ----
  dump      <- LoadGradFieldDump(data_dump$folder, data_dump$dumpfile)
  dump$dump <- DumpMesh(list_mesh$mesh, dump$dump)
  dump$wall <- DumpWall(list_mesh$wall, dump$dump)
  #--- * Accleration Data                                           ----
  # Note that tangent direction based on spline calc
  LoadSeshBCEq(data_dump$seshpath, "MOD_ALPHA_X")
  dump$a    <- BC_mod_alpha_x(dump$time)
  dump$wall <- DumpAccel(dump$a, dump$wall)
  #--- * Pressure Data                                              ----
  dump$wall <- DumpPres(dump$wall, interp = FALSE)
  #--- * Vorticity Data                                             ----
  dump$wall <- DumpVortOnly(dump$wall, dump$kinvis)
  #--- > Dump Calc Output                                           ----
  data_plot <- bind_rows(dump[c("time", "kinvis", "a")])
  data_plot <- cbind(data_dump, data_plot)
  list_dump <- c(
    data_plot = list(data_plot), 
    dump = dump[c("wall")])
  names(list_dump) <- c("data_plot", "wall")
  #--- > Produce Plots if Required                                  ----
  if (FALSE) {
    ggplot(dump$dump, aes(x, y, colour = u)) +
      geom_point() +
      scale_colour_gradientn(colours=rev(spectralpalette(10))) +
      coord_fixed(xlim=c(-0.4, 0.6), ylim=c(-0.2, 0.2))
    
    ggplot(dump$dump, aes(x, y, colour = o)) +
      geom_point() +
      scale_colour_gradientn(colours=rev(spectralpalette(10))) +
      coord_fixed(xlim=c(-0.4, 0.6), ylim=c(-0.2, 0.2))
    
    redu <- filter(dump$dump, x >= -0.5, x <= 0.8, y >= -0.2, y <= 0.2) %>% 
      select(x, y, u, v, p, o) %>% 
      unique(.)
    grid <- interp(x = round(redu$x, 10),
                   y = round(redu$y, 10),
                   z = round(redu$o, 10),
                   nx = 800,
                   ny = 600,
                   duplicate = "strip",
                   linear = TRUE,
                   extrap = FALSE)
    grid <- interp2xyz(grid, data.frame = TRUE)
    
    ggplot()  +
      geom_raster(aes(x, y, fill = z), grid, interpolate = TRUE) +
      geom_polygon(aes(x, y), dump$wall, fill = "white", colour = NA) +
      geom_path(aes(x, y, colour = o), dump$wall, size = 0.2) +
      scale_fill_gradientn(colours=rev(spectralpalette(10)),
                           limits = c(min(redu$o), max(redu$o)),
                           na.value = "white") +
      scale_colour_gradientn(colours=rev(spectralpalette(10)),
                           limits = c(min(redu$o), max(redu$o)),
                           na.value = "white",
                           guide = "none") +
      coord_fixed(xlim=c(-0.5, 0.8), ylim=c(-0.2, 0.2), expand = FALSE)
    
    
    
    x = c(-0.4, 0)
    
    velvec <- function(t, x) {
      u <- interp(
        x = round(redu$x, 12),
        y = round(redu$y, 12),
        z = round(redu$u, 12),
        xo = x[1],
        yo = x[2],
        duplicate = "strip",
        linear = TRUE)$z
      v <- interp(
        x = round(redu$x, 12),
        y = round(redu$y, 12),
        z = round(redu$v, 12),
        xo = x[1],
        yo = x[2],
        duplicate = "strip",
        linear = TRUE)$z
      return(c(u, v))
    }
    
    
    
  }
  
  rm(data_dump, data_plot, dump)
  # output
  return(list_dump)
}
