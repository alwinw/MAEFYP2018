#============================#
# NACA0012 Example File
# Alwin Wang
#----------------------------#

#--- Set Up ----
# Use rstudioapi to get saved location of this file; use str to print structures
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
# Source Required Scripts
source("src_library-manager.R")                                 # Call libraries and install missing ones
source("src_numerical-methods.R")                               # Load custom numerical methods
source("src_load-files.R")                                      # Load data
source("src_airfoil-analysis.R")                                # Airfoil files
source("src_vorticity-generation.R")                            # Vorticity Generation
source("src_plot-output.R")                                     # Plots to output and save
# Custom ggplot2 setup
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

# Path to example files for NACA0012
path = "../src-example/NACA0012/results/"
sesh = "RE-10000-sine-0.001-2000"
bndry <- LoadBndry(paste0(path, "bndry_prf"))
wall  <- LoadWallmsh(paste0(path, sesh))

# Plot airfoil normals
wallplot <- select(wall, -area) %>%
  mutate(nlen = sqrt(nx^2 + ny^2))
wallplot <- rbind(wallplot[,1:2], wallplot[,1:2] - wallplot[,3:4])
wallplot$wnum <- 1:nrow(wall)
ggplot(wallplot, aes(x, y, group = wnum)) +
  geom_line()
