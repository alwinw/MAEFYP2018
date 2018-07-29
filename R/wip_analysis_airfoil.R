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
saveplot   = "../src-example/NACA0012/"                         # Pass in later to plot_data
airfoil    = "NACA0012-AoA04"
folderpath = "../src-example/NACA0012/results/"
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

#--- * Session Data                                               ----
session   <- LoadSeshFileKeywords(data_mesh$seshpath)
long$sesh <- LongSesh(session)
plot_sesh <- ggplot(long$sesh, 
       aes(x, y, group=enum, colour=enum)) +
  geom_polygon(fill=NA) +
  geom_text(aes(elabx, elaby, label=enum, size=area),
            data = filter(long$sesh, ncorner=="n1")) +
  scale_color_gradientn(colours=spectralpalette(100))
# Plots
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
