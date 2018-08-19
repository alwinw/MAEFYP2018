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
# Plots
auxplot = 0

#--- Airfoil Calculation                                          ----
# Determine things like spline distance once per unique airfoil i.e. boundary profile and wall output
#--- * Boundary Data                                              ----
# Not sure this is actually used anywhere...
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
if (auxplot > 1) {
  # long_wallplot make normals etc
  ggplot(long_wall, aes(x, y, colour = theta, shape = up)) + 
    geom_point() + geom_path() + coord_fixed()
}
long_wall <- AirfoilSpline(long_wall)
long_wall <- AirfoilNorm  (long_wall)
if (auxplot > 0) {
  long_wallplot1        <- rbind(long_wall[,1:2], long_wall[,1:2] - long_wall[,3:4])
  long_wallplot1$wnum   <- 1:nrow(long_wall)
  long_wallplot1$type   <- "WallGrad"
  long_wallplot2        <- rbind(long_wall[,1:2], long_wall[,1:2] - long_wall[,16:17])
  long_wallplot2$vecerr <- long_wall$vecerr
  long_wallplot2$type   <- "Spline"
  ggplot(mapping = aes(x, y)) +
    geom_line(aes(group = wnum), linetype = "dashed", alpha = 0.5, data = long_wallplot1) +
    geom_point(aes(colour = vecerr), data = long_wallplot2) +
    scale_colour_gradientn(colours=spectralpalette(10))
  rm(long_wallplot1, long_wallplot2)
}
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
long$mesh <- LongMesh(long$mesh, long$sesh)
if (auxplot > 1) {
  ggplot(long$mesh,
         aes(x, y, group=enum, colour=enum)) +
    geom_point(shape='o', alpha = 0.2) +
    geom_polygon(fill=NA,
                 data = long$mesh %>% filter(node) %>% arrange(ncorner)) +
    geom_text(aes(elabx, elaby, label=enum, size=area),
              data = filter(long$mesh, ncorner=="n1")) +
    scale_color_gradientn(colours=spectralpalette(100)) +
    scale_size(guide="none", range=c(1*0.3, 6*0.8)) +
    coord_fixed()
}
#--- * Wall Data                                                  ----
long$wall <- list_airfoil$long_wall
long$wall <- LongWall(long$wall, long$mesh)
if (auxplot > 0) {
  long_wallplot <- long$wall %>%
    mutate(xd = x - nxG*aveh, yd = y - nyG*aveh) %>%
    select(x, y, xd, yd, enum, wnum)
  long_wallplot <- cbind(rbind(
    long_wallplot[,1:2], data.frame(x=long_wallplot$xd, y=long_wallplot$yd)),
    enum = long$wall$enum, wnum = long$wall$wnum)
  ggplot(long$wall,
         aes(x, y, group=enum, colour=enum)) +
    geom_point(aes(shape=node)) +
    geom_path() + 
    geom_path(aes(group=wnum), data=long_wallplot) + 
    scale_colour_gradientn(colours=spectralpalette(10)) +
    coord_fixed()
  rm(long_wallplot)
}
#--- * Local Data                                                 ----
long$mesh <- LocalMesh(long$mesh, long$wall)
if (auxplot > 0) {
  ggplot(long$mesh,
         aes(x, y, group=enum, colour=local)) +
    geom_point(shape='o', alpha = 0.2) +
    geom_polygon(fill=NA,
                 data = long$mesh %>% filter(node) %>% arrange(ncorner)) +
    geom_text(aes(elabx, elaby, label=enum, size=area),
              data = filter(long$mesh, ncorner=="n1")) +
    scale_color_gradientn(colours=rev(spectralpalette(10))) +
    scale_size(guide="none", range=c(1*0.3, 6*0.8)) +
    coord_fixed()
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
if (auxplot > 1) {
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
if (auxplot > 1) {
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
if (auxplot > 0) {
  ggplot(dump$wall, aes(x, y, colour = as)) +
    geom_point() +
    scale_colour_gradientn(colours=rev(spectralpalette(10))) +
    coord_fixed(xlim=c(-0.4, 0.6), ylim=c(-0.2, 0.2))
}
#--- * Pressure Data                                              ----
dump$wall <- DumpPres(dump$wall)
if (auxplot > 0) {
  ggplot(dump$wall, aes(x, y, colour = dpdsG)) +
    geom_point() +
    scale_colour_gradientn(colours=rev(spectralpalette(10))) +
    coord_fixed(xlim=c(-0.4, 0.6), ylim=c(-0.2, 0.2))
}
#--- * Vorticity Data                                             ----
order = 3 # Ideally, scale should be order/N_P
dump$offs <- AirfoilOffset(dump$wall, nsteps = order, scale = (order)/data_dump$tokenvalue)
if (auxplot > 1) {
  ggplot(dump$offs, aes(x, y, colour = norm, group = onum)) +
    geom_point() +
    geom_line() +
    coord_fixed()
}
dump$offs <- DumpVortInterp(dump$offs, dump$dump, 
                            linear = TRUE, extrap = FALSE, round = NULL)
dump$offs <- FiniteDiff(dump$offs, "o", order = order)
dump$offs <- DumpVortGrad(dump$offs)
if (auxplot > 0) {
  dodzfit <- lm(dodzS~dodzG, dump$offs)
  ggplot(dump$offs, aes(dodzG, dodzS)) + geom_point() + coord_fixed()+
    geom_abline(slope = 1, intercept = 0) + 
    stat_smooth(method = "lm", col = "blue") +
    labs(title = paste("Adj R2 = ",signif(summary(dodzfit)$adj.r.squared, 5),
                       "Intercept =",signif(dodzfit$coef[[1]],5 ),
                       " Slope =",signif(dodzfit$coef[[2]], 5),
                       " P =",signif(summary(dodzfit)$coef[2,4], 5)))
  rm(dodzfit)
}
dump$wall <- DumpVortJoin(dump$wall, dump$offs)
#--- > Dump Calc Output                                           ----
data_plot <- bind_rows(dump[c("time", "kinvis", "a")])
data_plot <- cbind(data_dump, data_plot)
list_dump <- c(
  data_plot = list(data_plot), 
  dump = dump[c("wall", "offs", "dump")])
names(list_dump) <- c("data_plot", "wall", "offs", "dump")
rm(data_dump, data_plot, dump, order)  

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


plot_nstheme <- ggplot(plot_wall, aes(s)) + 
  geom_vline(xintercept = as.numeric(plot_setup$vlines),        # Vertical lines for LE, TE, LE
             colour = "grey", linetype = "dashed") +
  geom_label(aes(x, y, label = labels), plot_setup$surf,        # Surface labels
             colour = "grey") +
  xlab("s") + 
  scale_x_continuous(breaks = plot_setup$xbreaks, 
                     labels = function(x) sprintf("%.2f", x)) +
  ylab(NULL) + ylim(c(-40, 30)) +
  scale_color_manual(
    name = "Legend",
    values = c("dp/ds" = "red", "dV/dt" = "blue", "LHS" = "purple", "RHS" = "purple"),
    labels = c(
      expression(-frac(1, rho)~frac(partialdiff*p, partialdiff*s)), 
      expression(-frac(partialdiff*V, partialdiff*t)), 
      expression(-
        bgroup("(",
               frac(1, rho)~frac(partialdiff*p, partialdiff*s) +
                 frac(partialdiff*V, partialdiff*t), ")")),
      expression(-nu~frac(partialdiff*omega, partialdiff*z))),
    guide = guide_legend(
      override.aes = list(
        linetype = c(rep("solid", 2), "dashed", "blank"),
        shape = c(rep(NA, 3), "O"),
        alpha = rep(1, 4)))) +
  theme(legend.key.size = unit(2.25, "lines"),
        legend.text.align = 0.5,
        legend.direction = "vertical", 
        legend.position = "right",
        legend.background = element_rect(colour = "black", size = 0.3))

plot_nsS <- plot_nstheme +
  geom_path(aes(y = -as, colour = "dV/dt")) +                   # Acceleration terms
  geom_path(aes(y = +dpdsS, colour = "dp/ds")) +                # Pressure field
  geom_path(aes(y = - as + dpdsS, colour = "LHS"),              # LHS, acceleration + pressure
            linetype = "dashed") +
  geom_point(aes(s, -dodzS*plot_data$kinvis, colour = "RHS"),    # RHS, v * dw/dz
             plot_offs, shape = "o", alpha = 0.3) +
  ggtitle(plot_setup$title)

plot_nsG <- plot_nstheme +
  geom_path(aes(y = -as, colour = "dV/dt")) +                   # Acceleration terms
  geom_path(aes(y = + dpdsG, colour = "dp/ds")) +               # Pressure field
  geom_path(aes(y = - as + dpdsG, colour = "LHS"),              # LHS, acceleration + pressure
            linetype = "dashed") +
  geom_point(aes(s, -dodzG*plot_data$kinvis, colour = "RHS"),    # RHS, v * dw/dz
             plot_offs, shape = "o", alpha = 0.3) +
  ggtitle(plot_setup$title)

print(plot_nsS)
print(plot_nsG)


#--- Look at integrals ----
ggplot(plot$dump %>% filter(enum == 100), 
       aes(x, y)) +
  geom_text(aes(label = jnum)) +
  geom_label(aes(label = ncorner)) +
  geom_point(aes(colour = ifelse(node, "node", "poly")), shape = "O") +
  coord_fixed()

N_P = 5

cont <- c(1:N_P, (1:N_P)*N_P, (N_P^2:1)[1:N_P], (N_P:1)*N_P - (N_P - 1))
cont <- c(unique(cont), cont[1])

inte_test <- plot$dump %>% filter(enum == 100)

