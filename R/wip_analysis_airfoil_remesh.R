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
adjm$wall$inpx <- cubicspline(adjm$uniw$s, adjm$uniw$x, xi = adjm$wall$ints)
adjm$wall$inpy <- cubicspline(adjm$uniw$s, adjm$uniw$y, xi = adjm$wall$ints)

adjm$edge <- read.table(paste0(saveplot, "remeshedge.dat"), header = TRUE)
temp <- left_join(
  rename(adjm$edge, nnum = nnum1), 
  unique(select(long$sesh, x, y, nnum)), 
  by = "nnum") %>% 
  rename(x1 = x, y1 = y)
adjm$edge <- cbind(adjm$edge, select(temp, x1, y1))
temp <- left_join(
  rename(adjm$edge, nnum = nnum2), 
  unique(select(long$sesh, x, y, nnum)), 
  by = "nnum") %>% 
  rename(x2 = x, y2 = y)
adjm$edge <- cbind(adjm$edge, select(temp, x2, y2))

adjm$edge <- adjm$edge %>% 
  mutate(intx = (x2 - x1)*spls + x1,
         inty = (y2 - y1)*spls + y1)

write.table(select(adjm$edge, spnum, intx, inty), 
            paste0(saveplot, "remeshedgeout.dat"),
            row.names = FALSE)

adjm$shifted <- adjm$wall %>% 
  filter(!is.na(dels)) %>% 
  arrange(anum) %>% 
  select(inpx, inpy)
adjm$shifted$newnum <- 
  c(986, 978, 977, 969, 968, 960, 959, 951, 950, 942,
    987, 995, 996, 1003, 1004, 1011, 1026, 1033, 1034, 1041, 1042, 1050)

adjm$shifted <- left_join(adjm$edge, adjm$shifted, by = "newnum")
adjm$shifted$shift <- 
  ifelse(!is.na(lag(adjm$shifted$inpx)) & is.na(adjm$shifted$inpx), "lag", "no")
adjm$shifted$shift <- 
  ifelse(!is.na(lead(adjm$shifted$inpx)) & is.na(adjm$shifted$inpx), "lead", adjm$shifted$shift)

adjm$shifted$inpx <- 
  ifelse(adjm$shifted$shift == "lag", adjm$shifted$intx + lag(adjm$shifted$inpx - adjm$shifted$intx, 1), adjm$shifted$inpx)
adjm$shifted$inpy <- 
  ifelse(adjm$shifted$shift == "lag", adjm$shifted$inty + lag(adjm$shifted$inpy - adjm$shifted$inty, 1), adjm$shifted$inpy)
adjm$shifted$inpx <- 
  ifelse(adjm$shifted$shift == "lead", adjm$shifted$intx + lead(adjm$shifted$inpx - adjm$shifted$intx, 1), adjm$shifted$inpx)
adjm$shifted$inpy <- 
  ifelse(adjm$shifted$shift == "lead", adjm$shifted$inty + lead(adjm$shifted$inpy - adjm$shifted$inty, 1), adjm$shifted$inpy)

adjm$shifted$shift <- 
  ifelse(!is.na(lag(adjm$shifted$inpx)) & is.na(adjm$shifted$inpx), "lag", "no")
adjm$shifted$shift <- 
  ifelse(!is.na(lead(adjm$shifted$inpx)) & is.na(adjm$shifted$inpx), "lead", adjm$shifted$shift)
adjm$shifted$inpx <- 
  ifelse(adjm$shifted$shift == "lag", adjm$shifted$intx + lag(adjm$shifted$inpx - adjm$shifted$intx, 2), adjm$shifted$inpx)
adjm$shifted$inpy <- 
  ifelse(adjm$shifted$shift == "lag", adjm$shifted$inty + lag(adjm$shifted$inpy - adjm$shifted$inty, 2), adjm$shifted$inpy)
adjm$shifted$inpx <- 
  ifelse(adjm$shifted$shift == "lead", adjm$shifted$intx + lead(adjm$shifted$inpx - adjm$shifted$intx, 2), adjm$shifted$inpx)
adjm$shifted$inpy <- 
  ifelse(adjm$shifted$shift == "lead", adjm$shifted$inty + lead(adjm$shifted$inpy - adjm$shifted$inty, 2), adjm$shifted$inpy)

adjm$shifted$shift <- 
  ifelse(!is.na(lag(adjm$shifted$inpx)) & is.na(adjm$shifted$inpx), "lag", "no")
adjm$shifted$shift <- 
  ifelse(!is.na(lead(adjm$shifted$inpx)) & is.na(adjm$shifted$inpx), "lead", adjm$shifted$shift)
adjm$shifted$inpx <- 
  ifelse(adjm$shifted$shift == "lag", adjm$shifted$intx + lag(adjm$shifted$inpx - adjm$shifted$intx, 3), adjm$shifted$inpx)
adjm$shifted$inpy <- 
  ifelse(adjm$shifted$shift == "lag", adjm$shifted$inty + lag(adjm$shifted$inpy - adjm$shifted$inty, 3), adjm$shifted$inpy)
adjm$shifted$inpx <- 
  ifelse(adjm$shifted$shift == "lead", adjm$shifted$intx + lead(adjm$shifted$inpx - adjm$shifted$intx, 3), adjm$shifted$inpx)
adjm$shifted$inpy <- 
  ifelse(adjm$shifted$shift == "lead", adjm$shifted$inty + lead(adjm$shifted$inpy - adjm$shifted$inty, 3), adjm$shifted$inpy)

adjm$shifted$inpx <- 
  ifelse(is.na(adjm$shifted$inpx), adjm$shifted$intx, adjm$shifted$inpx)
adjm$shifted$inpy <- 
  ifelse(is.na(adjm$shifted$inpy), adjm$shifted$inty, adjm$shifted$inpy)

write.table(select(adjm$shifted, newnum, inpx, inpy), 
            paste0(saveplot, "remeshedgeout.dat"),
            row.names = FALSE)

cat(
  sprintf("%5d %20.10E%20.10E%20.10E\n", adjm$shifted$newnum, adjm$shifted$inpx, adjm$shifted$inpy, rep(0, nrow(adjm$shifted))),
  file = paste0(saveplot, "remeshedgecat.dat"))


if (FALSE) {
  ggplot(long$sesh %>% filter(enum %in% adjm$adjm$enum), 
         aes(x, y, group = enum)) +
    geom_polygon(colour = "black", fill = NA)
  
  ggplot(long$wall %>% filter(enum %in% adjm$adjm$enum, node), 
         aes(x, y, group = enum)) +
    geom_point(colour = "black", shape = "O") +
    geom_point(aes(intx, inty), colour = "red",
               data = adjm$wall)
  
  ggplot(adjm$edge) +
    geom_point(aes(x1, y1), colour = "red", shape = "O") +
    geom_point(aes(x2, y2), colour = "red", shape = "O") +
    geom_point(aes(intx, inty), colour = "blue", shape = "O") +
    # geom_path(aes(intx, inty, group = sgrp), colour = "blue") +
    geom_text(aes(intx, inty, label = newnum), size = 3)
    geom_point(aes(inpx, inpy), colour = "green",
              data = adjm$shifted)
}













#--- * Required Paths                                             ----
airfoil    = "NACA0012_remesh"
folderpath = "../src-example/NACA0012_remesh/"
seshpath   = "NACA0012_remesh"
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

long$mesh <- LoadMesh(data_mesh$seshpath)
ggplot(long$mesh, aes(x, y, group = enum)) +
  geom_point() +
  # coord_cartesian(xlim = c(-0.45, -0.2), ylim = c(-0.1, 0.1))
  coord_cartesian(xlim = c(0.45, 0.7), ylim = c(-0.15, 0.1))
