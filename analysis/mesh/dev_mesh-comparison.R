#============================#
# Comparison of Meshes
# Alwin Wang
#----------------------------#

#--- Set Up                                                       ----
# Use rstudioapi to get saved location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
#--- * Scripts                                                    ----
# Source Required Scripts
srcpath = "../../R/"
source(paste0(srcpath, "src_library-manager.R"))                # Call libraries and install missing ones
# source(paste0(srcpath, "src_batch-functions.R"))                # Batch functions used
source(paste0(srcpath, "src_helper-functions.R"))               # Smaller functions used
# Additional scripts here
# ggplot2 setup (consider moving to a separate script)
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 

#--- * Function to Load Session Files                             ----
LoadSeshMesh <- function(session) {
  mesh <- LoadSeshFileKeywords(session)
  wall <- data.frame(
    enum = mesh$curves$element,
    node = TRUE)
  mesh <- LongSesh(mesh)
  mesh$node = TRUE
  mesh <- LocalMesh(mesh, wall)
  
  ggplot(mesh, aes(x, y, group = enum)) + 
    geom_polygon(aes(fill = local), colour = "white") +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(
      colours=rev(spectralpalette(max(mesh$local))))
  
  return(mesh)
}

naca0012  <- LoadSeshMesh("naca0012" )
naca0012r <- LoadSeshMesh("naca0012r")

plot_mesh <- ggplot() +
  geom_polygon(aes(x, y, group = enum, fill = local), 
               naca0012, colour = "white") +
  geom_polygon(aes(x, y, group = enum), 
               naca0012r, colour = "black", fill = NA) +
  geom_polygon(aes(x, y, group = enum), 
               naca0012, colour = "white", fill = NA) +
  geom_text(aes(elabx, elaby, label = enum, size = area),
            filter(naca0012, ncorner=="n1"), colour = "white") +
  geom_label(aes(x = x, y = y, label = nnum, size = area*0.2),
             naca0012r %>% ungroup() %>% group_by(nnum) %>% top_n(-1, area),
             label.padding = unit(0.1, "lines")) +
  scale_fill_gradientn(colours = brewer.pal(9, "Purples"),
                       limits = c(0, 9), guide = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot_mesh + coord_fixed(                                        # Trailing edge
  xlim = c(0.40, 0.70), ylim = c(-0.16, 0.06), expand = FALSE) +
  scale_size(guide="none", range=c(1*1, 6*8))
ggsave(paste0("mesh-te.png"),
       scale = 2, width = 10, height = 8, units = "cm", dpi = 300)
plot_mesh + coord_fixed(                                        # Leading edge
  xlim = c(-0.5, -0.2), ylim = c(-0.11, 0.13), expand = FALSE) +
  scale_size(guide="none", range=c(1*1, 6*8))
ggsave(paste0("mesh-le.png"),
       scale = 2, width = 10, height = 8, units = "cm", dpi = 300)
plot_mesh + coord_fixed(                                        # Airfoil
  xlim = c(-0.475, 0.675), ylim = c(-0.4, 0.4), expand = FALSE) +
  scale_size(guide="none", range=c(1*0.4, 6*4))


