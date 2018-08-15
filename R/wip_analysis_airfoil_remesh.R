#============================#
# Analysis of NACA0012 Example
# Alwin Wang
#----------------------------#

#--- Set Up                                                       ----
# Use rstudioapi to get saved location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
#--- * Scripts                                                    ----
# Source Required Scripts
source("src_library-manager.R")                                 # Call libraries and install missing ones
source("src_helper-functions.R")                                 # Smaller functions used
# Additional scripts here
# ggplot2 setup (consider moving to a separate script)
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 
#--- * Required Paths                                             ----
saveplot   = "../src-example/NACA0012_remesh/"                  # Pass in later to plot_data
airfoil    = "NACA0012-AoA04"
folderpath = "../src-example/NACA0012_remesh/results/"
seshpath   = "RE-10000-sine-0.001-2000"
dumppath   = "RE-10000-sine-0.001-2000-11.dump"
# Required input dataframes
data_airfoil <- data.frame(
  airfoil  = airfoil, 
  seshname = seshpath,
  folder   = folderpath,
  seshpath = paste0(folderpath, seshpath),
  stringsAsFactors = FALSE)
data_mesh <- data.frame(
  data_airfoil,
  tokenword  = "N_P",
  tokenvalue = 8,
  ID         = "NACA0012-AoA04-N_P5",
  stringsAsFactors = FALSE)
data_dump <- data.frame(
  data_mesh,
  dumpfile = dumppath,
  dumppath = paste0(folderpath, dumppath),
  stringsAsFactors = FALSE)
rm(airfoil, folderpath, seshpath, dumppath)

#--- Airfoil Calculation                                          ----
#--- * Boundary Data                                               ----
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
rm(data_airfoil, bndry, long_wall)

#--- Session and Mesh Calculation                                 ----
#--- * Airfoil Data                                               ----
long <- list()
#--- * Session Data                                               ----
session   <- LoadSeshFileKeywords(data_mesh$seshpath)
long$sesh <- LongSesh(session)
# Plots
if (FALSE) {
  # Plot setup
  plot_sesh <- ggplot(long$sesh, 
                      aes(x, y, group=enum, colour=enum)) +
    geom_polygon(fill=NA) +
    geom_text(aes(elabx, elaby, label=enum, size=area),
              data = filter(long$sesh, ncorner=="n1")) +
    scale_color_gradientn(colours=spectralpalette(10))
  # Plot output
  plot_sesh + coord_fixed() +                                     # Whole airfoil
    scale_size(guide="none", range=c(1*0.3, 6*0.8))
  plot_sesh + coord_fixed(                                        # Leading edge
    xlim = c(-0.45, -0.25), ylim = c(-0.055, 0.085), expand = FALSE) +
    scale_size(guide="none", range=c(1*3, 6*10))
  plot_sesh + coord_fixed(                                        # Trailing edge
    xlim = c(0.5, 0.675), ylim = c(-0.12, 0.024), expand = FALSE) +
    scale_size(guide="none", range=c(1*3, 6*10))
  
  plot_sesh + coord_fixed(                                        # Airfoil
    xlim = c(-0.45, 0.675), ylim = c(-0.3, 0.3), expand = FALSE) +
    scale_size(guide="none", range=c(1*0.5, 6*3))
}

#--- * Mesh Data                                                  ----
long$mesh <- LoadMesh(data_mesh$seshpath)
long$mesh <- LongMesh(long$mesh, long$sesh)
if (FALSE) {
  # Plot setup
  plot_mesh <- ggplot(long$mesh,
                      aes(x, y, group=enum, colour=enum)) +
    geom_point(shape='o', alpha = 0.2) +
    geom_polygon(fill=NA,
                 data = long$mesh %>% filter(node) %>% arrange(ncorner)) +
    geom_text(aes(elabx, elaby, label=enum, size=area),
              data = filter(long$mesh, ncorner=="n1")) +
    scale_color_gradientn(colours=spectralpalette(20))
  # Plot outputs
  plot_mesh + coord_fixed(                                        # Trailing edge
    xlim = c(0.4, 0.7), ylim = c(-0.15, 0.03), expand = FALSE) +
    scale_size(guide="none", range=c(1*3, 6*10))
}

#--- * Wall Data                                                  ----
long$wall <-  list_airfoil$long_wall
long$wall <- LongWall(long$wall, long$mesh)
#--- * Local Data                                                 ----
long$mesh <- LocalMesh(long$mesh, long$wall)
if (TRUE) {
  # Plot setup
  plot_mesh <- ggplot(long$mesh,
                      aes(x, y, group=enum, colour=local)) +
    geom_point(shape='o', alpha = 0.2) +
    geom_polygon(fill=NA,
                 data = long$mesh %>% filter(node) %>% arrange(ncorner)) +
    geom_text(aes(elabx, elaby, label=enum, size=area),
              data = filter(long$mesh, ncorner=="n1")) +
    scale_color_gradientn(colours=spectralpalette(20))
  # Plot outputs
  plot_mesh + coord_fixed(                                        # Trailing edge
    xlim = c(0.4, 0.8), ylim = c(-0.2, 0.045), expand = FALSE) +
    scale_size(guide="none", range=c(1*2, 6*8))
  plot_mesh + coord_fixed(                                        # Leading edge
    xlim = c(-0.5, -0.2), ylim = c(-0.10, 0.13), expand = FALSE) +
    scale_size(guide="none", range=c(1*2, 6*8))
  plot_mesh + coord_fixed(                                        # Airfoil
    xlim = c(-0.475, 0.675), ylim = c(-0.4, 0.4), expand = FALSE) +
    scale_size(guide="none", range=c(1*0.8, 6*4))
}
