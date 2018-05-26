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
# if (!file.info(savedata)$isdir) dir.create(savedata, recursive = TRUE)
# logfile = paste0(format(Sys.time(), "%Y-%m-%dT%H.%M.%S"), ".txt")

#--- List of Session FIles                                        ----
batchfolder = "../session-output"                               # Path to session files
batchdf <- ListSesh(batchfolder)                                # Dataframe of session files in folder

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
# Return the output as a list
airfoillistout <- list(
  airfoil = airfoilval, 
  bndry = bndry, 
  long_wall = long_wall)

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
# Offset
long$offset <- AirfoilOffset(long,                              # Create df of offset points from surface
                             totdist = 0.008, 
                             varh = TRUE, scale = 1)
long <- AirfoilOffsetEnum(long)                                 # Update enum values of the offset points
# Clean up and Return
meshlistout <- long[c("walldata", "threaddata", "offset")]

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
dumpval <- dumplist[[350]]
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
dump$offset = long$offset                                       # Initialise the offset for interpolation
dump <- DumpVortElements(dump, localval = 2, var = "t")         # Interpolate using each element

# Determine inflow velocity
# Determine Cp graph

# Derivatives
dump$offset <- FiniteDiff(dump$offset, "t_enum")                # Calculate derivative using each element
#--- Plot Output                                                  ----
PlotAccel(dump, dumpval,                                        # Plot acceleration curve
          save = FALSE)
PlotNS(dump, dumpval, long,                                     # Plot N-S for LHS vs RHS 
       save = FALSE)
PlotLE(dump, dumpval,                                           # Plot leading edge
       save = FALSE)