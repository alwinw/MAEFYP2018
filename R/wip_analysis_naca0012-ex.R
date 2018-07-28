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
#--- * Required Paths                                           ----
saveplot   = "../src-example/NACA0012/"
airfoil    = "NACA0012-AoA04"
folderpath = "../src-example/NACA0012/results/"
seshpath   = "RE-10000-sine-0.001-2000"
dumppath   = "RE-10000-sine-0.001-2000-03.dump"
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
  tokenvalue = 5,
  ID         = "NACA0012-AoA04-N_P5",
  stringsAsFactors = FALSE)
data_dump <- data.frame(
  data_mesh,
  dumpfile = dumppath,
  dumppath = paste0(folderpath, dumppath),
  stringsAsFactors = FALSE)
rm(airfoil, folderpath, seshpath, dumppath)
# Plots
auxplot = TRUE

#--- Airfoil Calculation                                          ----
# Determine things like spline distance once per unique airfoil i.e. boundary profile and wall output
#--- * Boundary Data                                              ----
bndry <- LoadBndry(data_airfoil$folder)
if (auxplot) {
  ggplot(bndry, aes(x, y)) + geom_point() + coord_fixed()
}
#--- * Wall Mesh Data                                             ----
wallmesh <- LoadWallGrad(data_airfoil$seshpath)
if (auxplot) {
  wallmeshplot      <- rbind(wallmesh[,1:2], wallmesh[,1:2] - wallmesh[,3:4])
  wallmeshplot$wnum <- 1:nrow(wallmesh)
  ggplot(wallmeshplot, aes(x, y, group = wnum)) + geom_line()
  rm(wallmeshplot)
}
long_wall <- AirfoilLongWall(wallmesh)
if (auxplot) {
  ggplot(long_wall, aes(x, y, colour = theta, shape = up)) + 
    geom_point() + geom_path() + coord_fixed()
}
long_wall <- AirfoilSpline(long_wall)
#--- > Airfoil Calc Output                                        ----
list_airfoil <- list(
  airfoil   = data_airfoil,
  bndry     = bndry,
  long_wall = long_wall)
rm(data_airfoil, bndry, long_wall)

#--- Session and Mesh Calculation                                 ----
#--- * Airfoil Data                                               ----
long <- list()
long$walldata = list_airfoil$long_wall
#--- * Session Data                                               ----
session       <- LoadSeshFileKeywords(data_mesh$seshpath)
long$seshdata <- 

if (auxplot) {
  ggplot(session)
}
