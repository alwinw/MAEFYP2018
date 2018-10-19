#============================#
# Analysis of NACA0012 Example
# Alwin Wang
#----------------------------#

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
saveplot   = "../src-example/Vincent/"                          # Pass in later to plot_data
airfoil    = "NACA0012-Vincent-2014"
folderpath = "../src-example/Vincent/results/"
seshpath   = "vincent"
dumppath   = "vincent-05.dump"
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

#--- > Sesh & Mesh Calc Output                                    ----
list_mesh <- list(
  wall = long$wall,
  mesh = long$mesh)
rm(data_mesh, long, session, wallmesh)

#--- Dump File Calculation                                        ----
#--- * Dump Data                                                  ----
dump      <- LoadGradFieldDump(data_dump$folder, data_dump$dumpfile)
dump$dump <- DumpMesh(list_mesh$mesh, dump$dump)
if (FALSE) {
  ggplot(dump$dump, aes(x, y, colour = u)) +
    geom_point() +
    scale_colour_gradientn(colours=rev(spectralpalette(10))) +
    coord_fixed(xlim=c(-0.4, 0.6), ylim=c(-0.2, 0.2))
  ggplot(dump$dump, aes(x, y, colour = p)) +
    geom_point() +
    scale_colour_gradientn(colours=rev(spectralpalette(10))) +
    coord_fixed(xlim=c(-0.4, 0.6), ylim=c(-0.2, 0.2))
  ggplot(dump$dump, aes(x, y, colour = o)) +
    geom_point() +
    scale_colour_gradientn(colours=rev(spectralpalette(10))) +
    coord_fixed(xlim=c(-0.4, 0.6), ylim=c(-0.2, 0.2))
}
dump$wall <- DumpWall(list_mesh$wall, dump$dump)
if (FALSE) {
  ggplot(dump$wall, aes(x, y, colour = dpdx)) +
    geom_point() +
    scale_colour_gradientn(colours=rev(spectralpalette(10))) +
    coord_fixed(xlim=c(-0.4, 0.6), ylim=c(-0.2, 0.2))
}
#--- * Accleration Data                                           ----
# Note that tangent direction based on spline calc
LoadSeshBCEq(data_dump$seshpath, "MOD_ALPHA_X")
dump$a    <- BC_mod_alpha_x(dump$time)
dump$wall <- DumpAccel(dump$a, dump$wall)
if (FALSE) {
  ggplot(dump$wall, aes(x, y, colour = as)) +
    geom_point() +
    scale_colour_gradientn(colours=rev(spectralpalette(10))) +
    coord_fixed(xlim=c(-0.4, 0.6), ylim=c(-0.2, 0.2))
}
#--- * Pressure Data                                              ----
dump$wall <- DumpPres(dump$wall)
if (FALSE) {
  ggplot(dump$wall, aes(x, y, colour = dpdsG)) +
    geom_point() +
    scale_colour_gradientn(colours=rev(spectralpalette(10))) +
    coord_fixed(xlim=c(-0.4, 0.6), ylim=c(-0.2, 0.2))
}
#--- * Vorticity Data                                             ----
dump$wall <- DumpVortOnly(dump$wall, dump$kinvis)

#--- > Dump Calc Output                                           ----
data_plot <- bind_rows(dump[c("time", "kinvis", "a")])
data_plot <- cbind(data_dump, data_plot)
list_dump <- c(
  data_plot = list(data_plot), 
  dump = dump[c("wall", "offs", "dump")])
names(list_dump) <- c("data_plot", "wall", "offs", "dump")
rm(data_dump, data_plot, dump)  

#--- Plot Outputs                                                 ----
#--- * Plot Setup                                                 ----
plot <- list()
plot$data      <- list_dump$data_plot
plot$data$save <- saveplot                                      # saveplot needs to be passed in!
plot$wall      <- list_dump$wall
plot$offs      <- list_dump$offs
plot$dump      <- list_dump$dump
plot$setup     <- PlotSetup(plot$wall, plot$data)
#--- * Plot Navier Stokes LHS vs RHS                          ----

plot_data <- plot$data
plot_wall <- plot$wall
plot_offs <- plot$offs
plot_dump <- plot$dump
plot_setup <- plot$setup

plot_title <- paste0(
  plot_data$airfoil, "\n",
  paste("Time:",                sprintf("%05.3f",  plot_data$time  )), "   ",
  paste("Acceleration:",        sprintf("%+07.4f", plot_data$a   )), "\n",
  paste("Kinematic Viscosity:", sprintf("%.4f",   plot_data$kinvis))
)


plot_nstheme <- ggplot(plot_wall, aes(s)) + 
  geom_vline(xintercept = as.numeric(plot_setup$vlines),        # Vertical lines for LE, TE, LE
             colour = "grey", linetype = "dashed") +
  # geom_label(aes(x, y, label = labels), plot_setup$surf,        # Surface labels
  geom_label(aes(x, y+20, label = labels), plot_setup$surf,        # Surface labels
             colour = "grey") +
  xlab("s") + 
  scale_x_continuous(breaks = plot_setup$xbreaks, 
                     labels = function(x) sprintf("%.2f", x)) +
  # ylab(NULL) + ylim(c(-40, 30)) +
  ylab(NULL) + ylim(c(-20, 20)) +
  scale_color_manual(
    name = "Legend",
    values = c("dp/ds" = "red", "dV/dt" = "blue", "LHS" = "purple", "RHS" = "purple"),
    labels = c(
      expression(-frac(1, rho)~frac(partialdiff*p, partialdiff*s)), 
      expression(-frac(partialdiff*V, partialdiff*t)), 
      expression(-bgroup("(",
                  frac(1, rho)~frac(partialdiff*p, partialdiff*s) +
                  frac(partialdiff*V, partialdiff*t), ")")),
      expression(-nu~frac(partialdiff*omega, partialdiff*z))),
    guide = guide_legend(
      override.aes = list(
        linetype = c(rep("solid", 2), "dashed", "blank"),
        # shape = c(rep(NA, 3), "O"),
        shape = c(rep(NA, 3), 20),
        alpha = rep(1, 4)))) +
  theme(legend.key.size = unit(2.25, "lines"),
        legend.text.align = 0.5,
        legend.direction = "vertical", 
        legend.position = "right",
        legend.background = element_rect(colour = "black", size = 0.3))

plot_nsG <- plot_nstheme +
  geom_path(aes(y = -as, colour = "dV/dt")) +                   # Acceleration terms
  geom_path(aes(y = + dpdsG, colour = "dp/ds")) +               # Pressure field
  geom_path(aes(y = - as + dpdsG, colour = "LHS"),              # LHS, acceleration + pressure
            linetype = "dashed") +
  geom_point(aes(s, -dodzG*plot_data$kinvis, colour = "RHS"),    # RHS, v * dw/dz
             plot_offs, shape = 20) +
             # plot_offs, shape = "o", alpha = 0.3) +
  ggtitle(plot_setup$title)

# print(plot_nsS)
print(plot_nsG)

ggsave("NS_Vincent_NP5.png", plot = plot_nsG, scale = 1.5,
       width = 15, height = 8, units = "cm", dpi = 300)
