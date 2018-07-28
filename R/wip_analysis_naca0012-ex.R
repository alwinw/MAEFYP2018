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
auxplot = 1

#--- Airfoil Calculation                                          ----
# Determine things like spline distance once per unique airfoil i.e. boundary profile and wall output
#--- * Boundary Data                                              ----
bndry <- LoadBndry(data_airfoil$folder)
if (auxplot > 0) {
  ggplot(bndry, aes(x, y)) + geom_point() + coord_fixed()
}
#--- * Wall Mesh Data                                             ----
wallmesh <- LoadWallGrad(data_airfoil$seshpath)
if (auxplot > 1) {
  wallmeshplot      <- rbind(wallmesh[,1:2], wallmesh[,1:2] - wallmesh[,3:4])
  wallmeshplot$wnum <- 1:nrow(wallmesh)
  ggplot(wallmeshplot, aes(x, y, group = wnum)) + geom_line()
  rm(wallmeshplot)
}
long_wall <- AirfoilLongWall(wallmesh)
if (auxplot > 0) {
  # long_wallplot make normals etc
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
long$wall= list_airfoil$long_wall
#--- * Session Data                                               ----
session   <- LoadSeshFileKeywords(data_mesh$seshpath)
long$sesh <- LongSesh(session)
if (auxplot > 1) {
  ggplot(long$sesh, 
         aes(x, y, group=enum, colour=enum)) +
    geom_polygon(fill=NA) +
    geom_text(aes(elabx, elaby, label=enum, size=area),
              data = filter(long$sesh, ncorner=="n1")) +
    scale_color_gradientn(colours=spectralpalette(100)) +
    scale_size(guide="none", range=c(1*0.3, 6*0.8)) +
    coord_fixed()
}
#--- * Mesh Data                                                  ----
long$mesh <- LoadMesh(data_mesh$seshpath)
long$mesh <- LongMesh(long$mesh)
if (auxplot > 1) {
  mid = floor(median(long$mesh$jnum))
  ggplot(long$mesh,
         aes(x, y, group=enum, colour=jnum)) +
    geom_point(shape = 'o', alpha = 0.5) +
    geom_text(aes(x, y, label=enum),
              data = filter(long$mesh, jnum==mid)) +
    coord_fixed()
}
#--- * Thread Data                                                ----
long$mesh <- LongJoin(long$mesh, long$sesh)
if (auxplot > 0) {
  ggplot(long$mesh,
         aes(x, y, group=enum)) +
    geom_polygon(fill=NA,
                 data = filter(long$mesh, !is.na(nnum)))
}

# Then determine for the wall what the equivalent mesh elements are