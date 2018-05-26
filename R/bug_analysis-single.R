#============================#
# Single Analysis
# Alwin Wang
#----------------------------#

# Fix derivative terms!

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
# Output Location
saveplot = "../plot-output"
# savedata = "Output_Data"
if (!dir.exists(saveplot)) dir.create(saveplot, recursive = TRUE)

#--- List of Session FIles                                        ----
batchfolder = "../session-output"                               # Path to session files
batchdf <- ListSesh(batchfolder)                                # Dataframe of session files in folder
# Clean up
rm(batchfolder)

#--- Airfoil Calculation                                          ----
# Iterate over list of airfoils list to load data
airfoillist <- split(batchdf, batchdf$airfoil)                  # Determine unique airfoil types
airfoillist <- lapply(airfoillist, function(x) x[1,])           # Take only the first row in each
# Function to load airfoil data
airfoilval <- airfoillist[[1]]                                  # HARD CODED FOR SINGLE ANALYSIS!!!!!!!!!!!!
# Boundary Data
bndrypath = paste0(airfoilval$folder, "bndry_prf")              # Path to bndry file
bndry <- LoadBndry(bndrypath)                                   # Read the bndry file
# Wall Mesh Data
wallmsh <- LoadWallmsh(airfoilval$seshpath)                     # Read the wallmesh file
long_wall <- AirfoilLongWall(wallmsh)                           # Airfoil data --> long_wall
long_wall <- AirfoilSpline(long_wall)                           # Determine spline distance

# Add a leading edge and trailing edge dataframe and s for upper LE --> upper TE and lower LE --> lower TE
# Maybe derotate the chord line to get horizontal location along the airfoil

# Return the output as a list
airfoillistout <- list(
  airfoil = airfoilval, 
  bndry = bndry, 
  long_wall = long_wall)
# Clean up
rm(airfoilval, bndry, long_wall, wallmsh)

#--- Session and Mesh Calculation                                 ----
# Load session info e.g. tokenwords = list("N_P", "N_Z")
tokenwords = list("N_P")                                        # Find unique bndry and knot N values
meshdf <- batchdf %>%
  rowwise() %>%
  do(data.frame(
    ., LoadSeshTokenWords(.$seshpath, tokenwords),
    stringsAsFactors = FALSE)) %>%
  mutate(ID = paste(
    airfoil, paste0(tokenword, tokenvalue), sep = "-"))  # This CANNOT handle multiple token words...
meshlist <- split(meshdf, meshdf$ID)                            # Unique airfoil and N_P
meshlist <- lapply(meshlist, function(x) x[1,])                 # Take only the first row in each
# Function to load mesh data
meshval <- meshlist[[2]]
long <- list()                                                  # Create a list of long format data
# Airfoil Data
airfoildata = airfoillistout                                    # HARD CODED FOR SINGLE ANALYSIS!!!!!!!!!!!!
long$walldata = airfoildata$long_wall                           # Airfoil --> long_walldata
# Session Data
keywords <- list(                                               # Keywords in session to read
  c("NODES", "nnum", "x", "y", "z"),
  c("ELEMENTS", "enum", "shapetag", "n1", "n2", "n3", "n4", "junk"),
  c("SURFACES", "snum", "element", "side", "bctag", "bc", "junk"),
  c("CURVES", "cnum", "element", "side", "curvetag", "curvedata", "junk"))
session <- LoadSeshFileKeywords(meshval$seshpath, keywords)     # Read keywords from session file
long$seshdata <- LongSession(session)                           # Session --> long_seshdata
# Mesh Data
mesh <- LoadMesh(meshval$seshpath)                              # Load mesh from mesh file
long$meshdata <- LongMesh(mesh)                                 # Mesh data --> long_meshdata
# Join Data
long$threaddata = long$meshdata                                 # Start with LARGEST data
long$threaddata = LongJoin(long$threaddata, long$seshdata)      # Left join with session data
long$threaddata <- LongAirfoil(long)                            # Special join for airfoil data
long$threaddata$wall = !is.na(long$threaddata$wnum)
long$threaddata$seshnode = !is.na(long$threaddata$nnum)
# Local Mesh
long <- LocalWallH(long)                                        # Average height of elements at wall
long$threaddata <-  LongJoin(                                   # Join with threaddata
  long$threaddata, select(long$wallaveh, enum, base, aveh, ar))
long$local <- LocalMesh(long)                                   # Determine how local the sessions elements are
long$threaddata <- LongJoin(
  long$threaddata, select(long$local, -nnum))
# OFFSET MOVED TO DUMP CALCULATION
# Clean up and Return
meshlistout <- long[c("walldata", "threaddata", "offset")]
rm(airfoildata, keywords, long, mesh, meshlist, meshval, session, tokenwords)

#--- Dump File List                                               ----
# Process batch list
dumpdf <- meshdf %>%                                            # Determine dump dataframe
  rowwise() %>% 
  do(data.frame(., dumpfile = ListDump(.$folder, .$seshname),
                stringsAsFactors = FALSE)) %>%
  mutate(dumppath = paste0(folder, dumpfile))
dumplist <- split(dumpdf, dumpdf$dumppath)                      # Create dump list
# Function to load and process dump files
# dumpval <- dumplist[[1202]]
dumpval <- dumplist[[391]]
# Long Mesh Data
long <- meshlistout                                             # HARD CODED FOR SINGLE ANALYSIS!!!!!!!!!!!!
# Dump Data
dump <- LoadDump(dumpval$folder, dumpval$dumpfile)              # Load dump file as list
dump$threaddata <- LongDump(dump, long)
# Acceleration
LoadSeshBCEqs(dumpval$seshpath, "MOD_ALPHA_X")                  # Load BC Equation from session file
dump$acceleration = BC_mod_alpha_x(dump$time)                   # Instantaenous acceleration
dump$accel = DumpAcceleration(dump, long)                       # Calculate acceleration components
dump$threaddata = LongJoin(dump$threaddata, dump$accel)         # Join with rest of data
# Pressure
dump$pres = DumpPressureStream(dump, long)                      # Calculate pressure gradient dp/ds
dump$threaddata = LongJoin(dump$threaddata, dump$pres)          # Join with rest of data
# Vorticity Interpolation
# Offset
long$offset <- AirfoilOffset(                                   # Create df of offset points from surface
  long, totdist = 0.008, varh = TRUE, scale = 0.25)
long <- AirfoilOffsetEnum(long)                                 # Update enum values of the offset points
dump$offset = long$offset                                       # Initialise the offset for interpolation
dump <- DumpVortElements(dump, localval = 3, var = "t")         # Interpolate using each element
# Derivatives
dump$offset <- FiniteDiff(dump$offset, "t_enum")                # Calculate derivative using each element

# Determine inflow velocity
# Determine Cp graph

surface <- dump$threaddata %>%                           # Initiate plot with threaddata
  filter(wall) %>% arrange(s) %>%
  select(-z, -elabx, -elaby, -area, -wall, -dxds, -dyds, -dydx,
         -dydxlen, -base, -aveh, -ar, -local) %>%
  mutate()
surface <- LongJoin(surface, dump$offset %>% select(s, t_enum_diff) %>% unique(.)) 

surface <- surface %>%
  mutate(LHS = - accels + dpds,
         RHS = t_enum_diff*dump$kinvis,
         error = (RHS - LHS)/RHS * 100)

velocity <- rbind(
  cbind(filter(dump$threaddata, x == min(dump$threaddata$x)), inflow = TRUE),
  cbind(filter(dump$threaddata, x == max(dump$threaddata$x)), inflow = FALSE))

ggplot(velocity, aes(u, y, group = inflow, colour = inflow)) +
  geom_point() +
  geom_path()
#+
  # facet_wrap(~ifelse(inflow, "inflow", "outflow"))

summary <- data.frame(
  x = c(2.5, 3),
  y = rep(-90, 2),
  labels = c(paste("median:", sprintf("%.2f%%", median(filter(surface, up)$error)), "\n",
                   "iqr:", sprintf("%.2f%%", IQR(filter(surface, up)$error))),
             paste("median:", sprintf("%.2f%%", median(filter(surface, !up)$error)), "\n",
                   "iqr:", sprintf("%.2f%%", IQR(filter(surface, !up)$error)))))

plot_vline <- data.frame(                                       # LE, TE, LE limits for vertical lines
  te_up = min(long$walldata$s),
  le = long$walldata[long$walldata$theta == pi,]$s[1],
  te_lo = max(long$walldata$s))
plot_surf <- data.frame(                                        # Labels for surfaces
  x = c(plot_vline$te_up,
        (plot_vline$te_up + plot_vline$le)/2,
        plot_vline$le,
        (plot_vline$le + plot_vline$te_lo)/2,
        plot_vline$te_lo),
  y = rep(-100, 5),
  labels = c("TE", "Upper Surface", "LE", "Lower Surface", "TE"))
plot_xbreaks <-                                                 # Breaks in x axis
  c(seq(plot_vline$te_up, plot_vline$le, length.out = 3)[1:2],
    seq(plot_vline$le, plot_vline$te_lo, length.out = 3))
plot_xlabs <- sprintf("%.2f", plot_xbreaks)
plot_xbreaks <- c(plot_xbreaks, 2.5, 3)
plot_xlabs <- c(plot_xlabs, "Uppper", "Lower")

ggplot(surface) +
  geom_vline(xintercept = as.numeric(plot_vline),               # Vertical lines for LE, TE, LE
             colour = "grey", linetype = "dashed") +
  geom_col(aes(s, error, fill = up)) +
  geom_boxplot(aes(as.numeric(up)*-0.5 + 3, error, colour = up),
               outlier.shape = "x") +
  geom_label(aes(x, y, label = labels), plot_surf,              # Surface labels
             colour = "grey40") +
  geom_label(aes(x, y, label = labels), summary,                # Summary labels
             colour = "grey40") +
  scale_x_continuous(breaks = plot_xbreaks, 
                     labels = plot_xlabs) +
  coord_cartesian(ylim = c(-100, 100)) +
  ylab(expression(frac(RHS - LHS, RHS)%*%100~"%")) +
  guides(fill = FALSE, colour = FALSE) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)))

ggplot(surface, aes(s, p, colour = up)) +
  geom_path()


#--- Plot N-S LHS vs RHS                                          ----
# Vertical Lines
plot_vline <- data.frame(                                       # LE, TE, LE limits for vertical lines
  te_up = min(long$walldata$s),
  le = long$walldata[long$walldata$theta == pi,]$s[1],
  te_lo = max(long$walldata$s))
plot_surf <- data.frame(                                        # Labels for surfaces
  x = c(plot_vline$te_up,
        (plot_vline$te_up + plot_vline$le)/2,
        plot_vline$le,
        (plot_vline$le + plot_vline$te_lo)/2,
        plot_vline$te_lo),
  y = rep(-20, 5),
  labels = c("TE", "Upper Surface", "LE", "Lower Surface", "TE"))
plot_xbreaks <-                                                 # Breaks in x axis
  c(seq(plot_vline$te_up, plot_vline$le, length.out = 3)[1:2],
    seq(plot_vline$le, plot_vline$te_lo, length.out = 3))
# Title
plot_title <- paste(
  dumpval$airfoil,
  paste("Kinematic Viscosity:", sprintf("%.4f", dump$kinvis)),
  paste("Time:", sprintf("%4.2f", dump$time)),
  paste("Acceleration:", sprintf("%+6.3f", dump$acceleration)),
  sep = "\n")
# Plot
plot_NS <- ggplot(dump$threaddata %>%                           # Initiate plot with threaddata
                    filter(wall) %>% arrange(s), aes(s)) +
  geom_vline(xintercept = as.numeric(plot_vline),               # Vertical lines for LE, TE, LE
             colour = "grey", linetype = "dashed") +
  geom_label(aes(x, y, label = labels), plot_surf,              # Surface labels
             colour = "grey") +
  geom_path(aes(y = -accels, colour = "dU/dt")) +               # Acceleration terms
  geom_path(aes(y = dpds, colour = "dp/ds")) +                  # Pressure field
  geom_path(aes(y = - accels + dpds, colour = "LHS"),           # LHS, acceleration + pressure
            linetype = "dashed") +
  geom_point(aes(s, t_enum_diff*dump$kinvis, colour = "RHS"),   # RHS, v * dw/dz
             dump$offset, shape = "O") +
  xlab("s") + 
  scale_x_continuous(breaks = plot_xbreaks, 
                     labels = function(x) sprintf("%.2f", x)) +
  ylab(NULL) + ylim(c(-20, 40)) +
  scale_color_manual(
    name = "Legend",
    values = c("dp/ds" = "red", "dU/dt" = "blue", "LHS" = "purple", "RHS" = "purple"),
    labels = c(
      expression(frac(1, rho)~frac(partialdiff*p, partialdiff*s)), 
      expression(frac(partialdiff*U, partialdiff*t)), 
      expression(
        bgroup("(",
               frac(1, rho)~frac(partialdiff*p, partialdiff*s) +
                 frac(partialdiff*U, partialdiff*t), ")")),
      expression(nu~frac(partialdiff*omega, partialdiff*z))),
    guide = guide_legend(
      override.aes = list(
        linetype = c(rep("solid", 2), "dashed", "blank"),
        shape = c(rep(NA, 3), "O")))) +
  theme(legend.key.size = unit(2.25, "lines"),
        legend.text.align = 0.5,
        legend.direction = "vertical", 
        legend.position = "right",
        legend.background = element_rect(colour = "black", size = 0.3)) +
  ggtitle(plot_title)
print(plot_NS)

#--- Plot Acceleration                                            ----
# ASSUME the time is 0 to 2 which may NOT always be the case
plot_aval <- data.frame(t = seq(0, 2, length.out = 500)) %>%
  mutate(a = BC_mod_alpha_x(t))
plot_accel <- ggplot(plot_aval, aes(t, a)) +
  geom_path(colour = "blue") +
  geom_point(aes(dump$time, dump$acceleration),
             colour = "blue")
print(plot_accel)


#--- Plot Leading Edge                                            ----
# Determine correct element for points for interpolation
n = 250
# plot_limits <- data.frame(
#   x = c(-0.55, 0.65),
#   y = c(-0.3, 0.3))
plot_limits <- data.frame(
  x = c(-0.455, -0.245),
  y = c(-0.07, 0.1))
pts_df = data.frame(
  x = rep(seq(plot_limits$x[1], plot_limits$x[2], length.out = n), each = n),
  y = rep(seq(plot_limits$y[1], plot_limits$y[2], length.out = n), times = n))
poly_df <- dump$threaddata %>% 
  filter(seshnode) %>%
  arrange(enum, ncorner) %>%
  select(x, y, enum)
pts_df <- PointinElement(pts_df, poly_df)

# Threaddata for enum
pts_list <- split(pts_df, pts_df$enum)
plot_full_df <- dump$threaddata %>%
  filter(enum %in% names(pts_list)) %>%
  select(x, y, t, enum, nnum, ncorner, local)
plot_full_list <- split(plot_full_df, plot_full_df$enum)

plot_le_interp <- lapply(names(pts_list), function(enum) {
  plot_le_interp <- as.data.frame(
    interpp(x = plot_full_list[[enum]]$x, y = plot_full_list[[enum]]$y, z = plot_full_list[[enum]]$t,
            xo = pts_list[[enum]]$x, yo = pts_list[[enum]]$y,
            linear = TRUE,
            duplicate = "strip"))
  return(plot_le_interp)
})

plot_le_interp <- bind_rows(plot_le_interp) %>%
  filter(!is.na(z))
colnames(plot_le_interp) <- c("x", "y", "t")

PlotScaleLimit <- function(t, max) ifelse(abs(t) > max, max*sign(t), t)
PlotContourBreak <- function(max, step) seq(-max, max, by = step)[seq(-max, max, by = step) != 0]

plot_le <- ggplot(plot_le_interp) +
  # Vorticity tile
  geom_tile(aes(x, y, fill = PlotScaleLimit(t, 200)), alpha = 0.8) +
  # Original Points
  # geom_point(aes(x, y, fill = PlotScaleLimit(t, 200)), 
             # shape = 23, colour = "white",
             # data = plot_full_df) +
  # Original Mesh
  geom_polygon(aes(x, y, group = enum), 
               fill = NA, colour = "grey90", alpha = 0.1, linetype = "1616",
               data = plot_full_df %>% filter(!is.na(nnum)) %>% arrange(enum, ncorner)) +
  # Vorticity contours
  # stat_contour(aes(x, y, z = t, fill = ..level..), geom = "polygon", alpha =1,
  # breaks = PlotContourBreak(200, 10)) +
  stat_contour(aes(x, y, z = t, colour = ..level..), alpha = 1,
               breaks = PlotContourBreak(200, 10)) +
  # Chosen offset points
  geom_point(aes(x, y, fill = PlotScaleLimit(t_enum, 200)),
               data = dump$offset,
               shape = 21) +
  # Airfoil
  geom_polygon(aes(x, y), fill = "white", colour = "black", alpha = 1,
               data = dump$accel) +
  coord_fixed(
    xlim = plot_limits$x,
    ylim = plot_limits$y,
    expand = FALSE) +
  scale_fill_gradientn(name = "vorticity\n", colours = spectralpalette(20), limits = c(-200, 200)) +
  scale_colour_gradientn(colours = spectralpalette(20), limits = c(-200, 200), guide = "none")
  # scale_colour_gradient(low = "grey99", high = "grey80", guide = "none")
print(plot_le)
