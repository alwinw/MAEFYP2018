#============================#
# NACA0012 Example File
# Alwin Wang
#----------------------------#

# May consider moving the smaller single functions to 3-5 larger functions that are run once

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
saveplot = "../src-example/NACA0012/"
path = "../src-example/NACA0012/results/"
sesh = "RE-10000-sine-0.001-2000"
dump = "RE-10000-sine-0.001-2000-03.dump"
airfoilval <- data.frame(
  folder = path, seshpath = paste0(path, sesh))
meshval <- data.frame(
  seshpath = paste0(path, sesh))
dumpval <- data.frame(
  folder = path, dumpfile = dump,
  seshpath = paste0(path, sesh))

#--- Airfoil Calculation                                          ----
#--- Boundary Data                                                ----
bndrypath = paste0(airfoilval$folder, "bndry_prf")              # Path to bndry file
bndry <- LoadBndry(bndrypath)                                   # Read the bndry file [x y]
# ggplot(bndry, aes(x, y)) + geom_point(shape = 'o') + coord_fixed()
# --- Wall Mesh Data                                               ----
wallmsh   <- LoadWallGrad(airfoilval$seshpath)                  # Read the wallgrad file [x y nxG nyG areaG]
# plot this showing the normals from the surface?
long_wall <- AirfoilLongWall(wallmsh)                           # LE->TE->LE [x y nxG nyG areaG wnum theta up wall]
long_wall <- AirfoilSpline(long_wall)                           # Spline dist [x y ... s dxds dydx dydx dydxlen]
check     <- AirfoilSplineCheck(long_wall)
# Return the output as a list
airfoildata <- list(
  airfoil = airfoilval, 
  bndry = bndry, 
  long_wall = long_wall)

#--- Session and Mesh Calculation                                 ----
long <- list()                                                  # Create a list of long format data
#--- Airfoil Data                                                 ----
long$walldata = airfoildata$long_wall                           # Airfoil --> long_walldata
#--- Session Data                                                 ----
keywords <- list(                                               # Keywords in session to read
  c("NODES", "nnum", "x", "y", "z"),
  c("ELEMENTS", "enum", "shapetag", "n1", "n2", "n3", "n4", "junk"),
  c("SURFACES", "snum", "element", "side", "bctag", "bc", "junk"),
  c("CURVES", "cnum", "element", "side", "curvetag", "curvedata", "junk"))
session <- LoadSeshFileKeywords(meshval$seshpath, keywords)     # Read keywords from session file
long$seshdata <- LongSession(session)                           # Session --> long_seshdata
#--- Mesh Data                                                    ----
mesh <- LoadMesh(meshval$seshpath)                              # Load mesh from mesh file
long$meshdata <- LongMesh(mesh)                                 # Mesh data --> long_meshdata
#--- Join Data                                                    ----
long$threaddata = long$meshdata                                 # Start with LARGEST data
long$threaddata = LongJoin(long$threaddata, long$seshdata)      # Left join with session data
long$threaddata <- LongAirfoil(long)                            # Special join for airfoil data
long$threaddata$wall = !is.na(long$threaddata$wnum)
long$threaddata$seshnode = !is.na(long$threaddata$nnum)
#--- Local Mesh                                                   ----
long <- LocalWallH(long)                                        # Average height of elements at wall
long$threaddata <-  LongJoin(                                   # Join with threaddata
  long$threaddata, select(long$wallaveh, enum, base, aveh, ar))
long$local <- LocalMesh(long)                                   # Determine how local the sessions elements are
long$threaddata <- LongJoin(
  long$threaddata, select(long$local, -nnum))
output <- long[c("walldata", "threaddata")]


#--- Dump File List                                               ----
#---  Long Mesh Data                                              ----
# long <- meshlistin[[dumpval$ID]]                                # Collect airfoil data
#--- Dump Data                                                    ----
dump <- LoadGradFieldDump(dumpval$folder, dumpval$dumpfile)     # Load dump file as list
dump$threaddata <- LongDump(dump, long)
#--- Acceleration                                                 ----
LoadSeshBCEqs(dumpval$seshpath, "MOD_ALPHA_X")                  # Load BC Equation from session file
dump$acceleration = BC_mod_alpha_x(dump$time)                   # Instantaenous acceleration
dump$accel = DumpAcceleration(dump, long)                       # Calculate acceleration components
dump$threaddata = LongJoin(dump$threaddata, dump$accel)         # Join with rest of data
#--- Pressure                                                     ----
dump$pres = DumpPressureStream(dump, long)                      # Calculate pressure gradient dp/ds
dump$threaddata = LongJoin(dump$threaddata, dump$pres)          # Join with rest of data
#--- Vorticity Interpolation                                      ----
# Offset
long$offset <- AirfoilOffset(                                   # Create df of offset points from surface
  long, totdist = 0.008, varh = TRUE, scale = 0.25)
long <- AirfoilOffsetEnum(long)                                 # Update enum values of the offset points
dump$offset = long$offset                                       # Initialise the offset for interpolation
# dump <- DumpVortTransformed(dump, localval = 2, var = "t")      # Interpolate based on stream, norm cs
dump <- DumpVortElements(dump, localval = 2, var = "t")         # Interpolate using each element

# Determine inflow velocity

# Derivatives
# dump$offset <- FiniteDiff(dump$offset, "t_trans")             # Calculate derivative using stream, norm cs
dump$offset <- FiniteDiff(dump$offset, "t_enum")                # Calculate derivative using each element
#--- Data Output                                                  ----
# Surface is the surface only data
surface <- dump$threaddata %>%
  filter(wall) %>% arrange(s) %>%                               # Only surface, only variables of interest
  select(-z, -elabx, -elaby, -area, -wall, -dxds, -dyds, -dydx,
         -dydxlen, -base, -aveh, -ar, -local)
surface <- LongJoin(surface,                                    # Join with dw/dz
                    dump$offset %>% select(s, t_enum_diff) %>% unique(.)) %>%
  mutate(LHS = - accels + dpds,
         RHS = t_enum_diff*dump$kinvis,
         error = (RHS - LHS)/RHS * 100)                         # Percentage error
# Velocity Profiles
velocity <- rbind(
  cbind(filter(dump$threaddata, x == min(dump$threaddata$x)), inflow = TRUE),
  cbind(filter(dump$threaddata, x == max(dump$threaddata$x)), inflow = FALSE))
dumpinfo <- data.frame(time = dump$time, kinvis = dump$kinvis, acceleration = dump$acceleration)
#--- Plot Output                                                  ----
# Plot acceleration curve
PlotAccel(dump, dumpval, save = TRUE, saveplot)
# Plot N-S for LHS vs RHS
PlotNS(dump, dumpval, long, save = TRUE, saveplot)
# Plot leading edge
PlotLE(dump, dumpval, save = TRUE, saveplot)
# Plot error
PlotError(dump, surface, dumpval, long, save = TRUE, saveplot)
# Return the output (turned off temporarily)
# return(list(dumpval = dumpval, surface = surface, velocity = velocity, dumpinfo = dumpinfo))


