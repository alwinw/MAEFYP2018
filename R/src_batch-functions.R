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
  #--- * Mass Data
  long$mass <- LoadMass(data_mesh$seshpath)
  long$mesh <- LongMass(long$mesh, long$mass)
  #--- * Wall Data                                                  ----
  long$wall <- list_airfoil$long_wall
  long$wall <- LongWall(long$wall, long$mesh)
  #--- * Local Data                                                 ----
  long$mesh <- LocalMesh(long$mesh, long$wall)
  # Determine enum regions as long$enum
  #--- > Sesh & Mesh Calc Output                                    ----
  list_mesh <- list(
    wall = long$wall,
    mesh = long$mesh)
  rm(data_mesh, long, session)
  # output
  return(list_mesh)
}

#--- Dump File Calculation                                        ----
BatchLoadDump <- function(data_dump, outp_mesh, plot = "none", outp = "wall",
                          srcpath = "", addscr = NULL) {
  source(paste0(srcpath, "src_library-manager.R"))                # Call libraries and install missing ones
  source(paste0(srcpath, "src_helper-functions.R"))               # Smaller functions used
  source(paste0(srcpath, "src_plotting-functions.R"))             # Functions for plotting
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
  #--- * Dump Calc Output                                           ----
  data_plot <- bind_rows(dump[c("time", "kinvis", "a")])
  data_plot <- cbind(data_dump, data_plot) %>% 
    mutate(plotname = paste0(
      airfoil, "-v", sprintf("%0.4f", kinvis), "-t", sprintf("%06.4f", time)))
  
  #--- * Integral Output                                            ----
  
  
  
  #--- * Create Output                                              ----
  if (outp == "wall") {
    list_dump <- c(
      data_plot = list(data_plot), 
      dump = dump[c("wall")])
    names(list_dump) <- c("data_plot", "wall")
  } else if (outp == "node") {
    dump$node <- filter(dump$dump, node)
    list_dump <- c(
      data_plot = list(data_plot), 
      dump = dump[c("wall")],
      node = dump[c("node")])
    names(list_dump) <- c("data_plot", "wall", "node")
  } else if (outp == "te") {
    xvec = c( 0.595,  0.615)
    yvec = c(-0.046, -0.030)
    xveci = xvec + c(-0.1*(xvec[2]-xvec[1]), 0.1*(xvec[2]-xvec[1]))
    yveci = yvec + c(-0.1*(yvec[2]-yvec[1]), 0.1*(yvec[2]-yvec[1]))
    dump$te <- filter(dump$dump, x >= xveci[1], x <= xveci[2], y >= yveci[1], y <= yveci[2]) %>% 
      select(x, y, u, v, p, o) %>% 
      unique(.)
    list_dump <- c(
      data_plot = list(data_plot), 
      dump = dump[c("wall")],
      node = dump[c("te")])
    names(list_dump) <- c("data_plot", "wall", "te")
  }
  
  #--- > Produce Plots if Required                                  ----
  # Last entry of plot needs to be saveplot
  if        ("none" %in% plot) {
  } else if ("NS" %in% plot) {
    
  } else {
    saveplot = plot[1]
    if        ("airfoil"  %in% plot) {
      # Plot the airfoil
      plot_af = PlotAirfoil(
        dump$dump, dump$wall, data_plot, scalearr = c(10, 4), saveplot) } 
    if ("TEstream" %in% plot) {
      # Plot the trailing edge
      plot_te = PlotTE(
        dump$dump, dump$wall, data_plot, scalearr = c(10, 4), saveplot) } 
    if ("LEstream" %in% plot) {
      
    }
  } 
  
  #--- > Run additional scripts if called                           ----
  if (!is.null(addscr))
    for (i in 1:length(addscr)) source(addscr[i])

  rm(data_dump, data_plot, dump)
  # output
  return(list_dump)
}
