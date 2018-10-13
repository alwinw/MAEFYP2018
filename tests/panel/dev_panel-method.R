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
panels <- VLPanel(coord)
ggplot(coord, aes(x, y)) + geom_path() + geom_point() + 
  coord_fixed() +
  geom_segment(aes(x, y, xend = x+nx, yend = y+ny), panels)
VLSol(coord, aoa = 5, U = 1)
actual <- data.frame(gamma = c(-0.8862, -0.9905, 0.4848, 1.233, 0.8862))

#--- * NACA2412 Airfoil                                           ----
# Source: Kuethe & Chow pg 121
coord <- data.frame(
  x = c(1.000, 0.933, 0.750, 0.500, 0.250, 0.067, 0.000, 0.067, 0.250, 0.500, 0.750, 0.933, 1.000),
  y = c(0.000,-0.005,-0.017,-0.033,-0.042,-0.033, 0.000, 0.045, 0.076, 0.072, 0.044, 0.013, 0.000) )
panels <- VLPanel(coord)
ggplot(coord, aes(x, y)) + geom_path() + geom_point() + 
  coord_fixed() +
  geom_segment(aes(x, y, xend = x+nx, yend = y+ny), panels)
VLSol(coord, aoa = 8, U = 1)

#--- * NACA0012 Airfoil                                           ----
# Source: Hugh Blackburn
coord <- data.frame(
  x = c(1.000000, 0.904508, 0.654508, 0.345491, 0.095491,
        0.000000, 0.095491, 0.345491, 0.654508, 0.904508,
        1.000000),
  y = c(0.000000,-0.013914,-0.040917,-0.059575,-0.046049,
        0.000000, 0.046049, 0.059575, 0.040917, 0.013914,
        0.000000) )
panels <- VLPanel(coord)
ggplot(coord, aes(x, y)) + geom_path() + geom_point() + 
  coord_fixed() +
  geom_segment(aes(x, y, xend = x+nx, yend = y+ny), panels)
VLSol(coord, aoa = 4, U = 1)