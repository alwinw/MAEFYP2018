#============================#
# Panel Method
# Alwin Wang
#----------------------------#

# Panel method: Vortex Panels of Linearly Varying Strength

#--- Set Up                                                       ----
# Use rstudioapi to get saved location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
#--- * Scripts                                                    ----
# Source Required Scripts
srcpath = "../../R/"
source(paste0(srcpath, "src_library-manager.R"))                # Call libraries and install missing ones
source(paste0("../../analysis/panel/", "src_panel-functions.R"))
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 


#--- Test Cases                                                   ----
#--- * Diamond Airfoil                                            ----
# Source: MAE4409 Lecture notes; Hugh Blackburn's C code
coord <- data.frame(
  x = c(1,  0.3,  0,  0.3, 1),
  y = c(0, -0.05, 0, 0.05, 0))
