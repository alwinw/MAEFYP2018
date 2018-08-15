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
#--- * Boundary Data                                              ----
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
if (FALSE) {
  # Plot setup
  plot_mesh <- ggplot(long$mesh,
                      aes(x, y, group=enum, colour=local)) +
    geom_point(shape='o', alpha = 0.2) +
    geom_polygon(fill=NA,
                 data = long$mesh %>% filter(node) %>% arrange(ncorner)) +
    geom_text(aes(elabx, elaby, label=enum, size=area),
              data = filter(long$mesh, ncorner=="n1")) +
    geom_label(aes(x = x, y = y, label = nnum, size = area*0.2), 
               label.padding = unit(0.1, "lines"),
              data = long$mesh %>% ungroup() %>% group_by(nnum) %>% top_n(-1, area)) +
    scale_color_gradientn(colours=spectralpalette(20),
                          guide = "none")
  # Plot outputs
  plot_mesh + coord_fixed(                                        # Trailing edge
    xlim = c(0.4, 0.8), ylim = c(-0.16, 0.07), expand = FALSE) +
    scale_size(guide="none", range=c(1*1, 6*8))
  plot_mesh + coord_fixed(                                        # Leading edge
    xlim = c(-0.5, -0.2), ylim = c(-0.10, 0.13), expand = FALSE) +
    scale_size(guide="none", range=c(1*1, 6*8))
  plot_mesh + coord_fixed(                                        # Airfoil
    xlim = c(-0.475, 0.675), ylim = c(-0.4, 0.4), expand = FALSE) +
    scale_size(guide="none", range=c(1*0.4, 6*4))
}

# Remesh Airfoil                                                  ----
# Note that this assumes the original airfoil is being used
# and that the direction of s is clockwise!

adjm      <- list()
# List of elements to be split
adjm$adjm <- data.frame(
  anum = seq(1, 60), 
  enum = c( 90,  57,  58,  59,  60, 381, 382, 383, 384, 361,  # LE local = 1
           116,  61,  62,  63,  64, 377, 378, 379, 380, 349,  # LE local = 2
           142, 149, 150, 151, 152, 373, 374, 375, 376, 337,  # LE local = 3
            66,  65, 144, 143, 704, 703, 336, 335, 372, 371,  # TE local = 1
            92,  91, 146, 145, 705, 702, 334, 333, 360, 359,  # TE local = 2
           118, 117, 148, 147, 706, 701, 332, 331, 348, 347), # TE local = 3
  local = rep(c(rep(1, 10), rep(2, 10), rep(3, 10)),  2),
  up    = rep(c(rep(TRUE, 5), rep(FALSE, 5)), 6),
  spls  = c(rep(c(0.4, rep(0.5, 8), 0.6), 3),
           rep(c(0.6, 0.5, rep(0.33, 2), rep(0.25, 2), rep(0.33, 2), 0.5, 0.4), 3)),
  spln  = rep(c(rep(0, 20), rep(0.5, 10)), 2)
)
# Double splits
adjm$doub <- adjm$adjm %>% 
  filter(spls %in% c(0.33, 0.25)) %>% 
  mutate(spls = spls + 0.33)
# Combine in double splits
adjm$adjm <- rbind(adjm$adjm, adjm$doub) %>% 
  arrange(anum, spls) %>% 
  mutate(anum = 1:n())
# Determine s value on surface
adjm$wall <- long$wall %>% 
  filter(enum %in% adjm$adjm$enum, node) %>% 
  select(x, y, s, enum, aveh) %>% 
  arrange(enum, s) %>% 
  group_by(enum) %>% 
  mutate(dels = lead(s) - s) %>% 
  ungroup()
adjm$wall <- left_join(adjm$wall, adjm$adjm, by = "enum") %>% 
  mutate(ints = s + spls*dels) %>% 
  group_by(anum) %>%
  mutate(ints = sum(ints, na.rm = TRUE)) %>% 
  ungroup() %>% 
  arrange(ints)
# Interpolate (x, y) on surface
adjm$uniw <- long$wall %>% 
  select(x, y, s)
adjm$uniw <- unique(adjm$uniw)
adjm$wall$intx <- cubicspline(adjm$uniw$s, adjm$uniw$x, xi = adjm$wall$ints)
adjm$wall$inty <- cubicspline(adjm$uniw$s, adjm$uniw$y, xi = adjm$wall$ints)


adjm$edge <- long$sesh %>%  
  filter(enum %in% adjm$adjm$enum) %>%  
  left_join(select(adjm$adjm, anum, enum), by = "enum") %>% 
  left_join(select(adjm$wall, x, y, local),
            by = c("x", "y")) %>% 
  unique(.) %>% 
  group_by(enum, anum) %>%  
  mutate(temp = sum(local, na.rm = TRUE)) %>%  
  mutate(local = ifelse(is.na(local) & temp == 2, 2, local)) %>% 
  select(-temp) %>% 
  filter(local %in% c(1, 2)) %>% 
  arrange(local, anum)


if (FALSE) {
  ggplot(long$sesh %>% filter(enum %in% adjm$adjm$enum), 
         aes(x, y, group = enum)) +
    geom_polygon(colour = "black", fill = NA)
  
  ggplot(long$wall %>% filter(enum %in% adjm$adjm$enum, node), 
         aes(x, y, group = enum)) +
    geom_point(colour = "black", shape = "O") +
    geom_point(aes(intx, inty), colour = "red",
               data = adjm$wall)
  
  ggplot(adjm$edge, aes(x, y, group = enum)) +
    geom_polygon(colour = "black", fill = NA) +
    geom_point(aes(intx, inty), colour = "red",
               data = adjm$wall)
}

