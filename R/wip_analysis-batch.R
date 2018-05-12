#============================#
# Batch Analysis
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
# Custom ggplot2 setup
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
#--- Load List of Session FIles                                   ----
batchfolder = "../session-files"                                # Path to session files
batchlist <- ListSesh(batchfolder)                              # List sessionfiles in folder

# Check if all file types exist, if not then call bash script

#--- Output Location ----
# saveplot = "Output_Plot"
# savedata = "Output_Data"
# if (!file.info(saveplot)$isdir) dir.create(saveplot, recursive = TRUE)
# if (!file.info(savedata)$isdir) dir.create(savedata, recursive = TRUE)
# logfile = paste0(format(Sys.time(), "%Y-%m-%dT%H.%M.%S"), ".txt")

#--- Airfoil Calculation                                          ----
airfoillist <- split(batchlist, batchlist$airfoil)              # Determine unique airfoil types
airfoillist <- lapply(airfoillist, function(x) x[1,])           # Take only the first row in each
# Function to load airfoil data
BatchLoadAirfoil <- function(airfoilval) {                      # airfoilval = airfoillist[[1]]
  source("src_load-files.R")                                      # Source required functions
  #--- Boundary Data                                                ----
  bndrypath = paste0(airfoilval$folder, "bndry_prf")              # Path to bndry file
  bndry <- LoadBndry(bndrypath)                                   # Read the bndry file
  #--- Wall Mesh Data                                               ----
  wallmsh <- LoadWallmsh(airfoilval$seshpath)                     # Read the wallmesh file
  long_wall <- AirfoilLongWall(wallmsh)                           # Airfoil data --> long_wall
  long_wall <- AirfoilSpline(long_wall)                           # Determine spline distance
  # unixy_wallmsh <-                                                # Unique (x,y) for left_join
    # long_wall[!duplicated(select(long_wall, x, y)),]
  # Return the output as a list
  return(list(
    airfoil = airfoilval, 
    bndry = bndry, 
    long_wall = long_wall))
    # , unixy_wallmsh = unixy_wallmsh))
}
# Iterate over list of airfoils list to load data
airfoillist <- lapply(airfoillist, BatchLoadAirfoil)

# cl <- makeCluster(detectCores())                                # Start the cluster
# clusterExport(cl, c("airfoillist"))                             # Export objects into the cluster
# # Function to be applied across the list
# #    pblappy here
# # Stop cluster
# stopCluster(cl)

#---- Dump File List                                              ----
dumplist <- batchlist %>%                                       # Determine dump list
  rowwise() %>% 
  do(data.frame(., dumpfile = ListDump(.$folder, .$seshname))) %>%
  mutate(dumppath = paste0(folder, "/", dumpfile))
threadlist <- split(dumplist, dumplist$dumppath)                # Create thread list

#--- Main Apply Function                                          ----
BatchThread <- function(threadval, airfoillist) {               # threadval = threadlist[[1]]
  # Assume:
  #  The order n x n may be different for each session file
  #  So, long$... must be generated for each session file
  #  even if the airfoil file is the same
  print(threadval$dumppath)
  long <- list()                                                  # Create a list of long format data
  #---- Airfoil Data                                                ----
  airfoildata = airfoillist[[threadval$airfoil]]                  # Collect airfoil data
  long$walldata = airfoildata$long_wall                           # Airfoil --> long_walldata
  #--- Session Data                                                 ----
  keywords <- list(                                               # Keywords in session to read
    c("NODES", "nnum", "x", "y", "z"),
    c("ELEMENTS", "enum", "shapetag", "n1", "n2", "n3", "n4", "junk"),
    c("SURFACES", "snum", "element", "side", "bctag", "bc", "junk"),
    c("CURVES", "cnum", "element", "side", "curvetag", "curvedata", "junk"))
  session <- LoadSeshFileKeywords(threadval$seshpath, keywords)   # Read keywords from session file
  LoadSeshBCEqs(threadval$seshpath, "MOD_ALPHA_X")                # Load BC Equation from session file
  long$seshdata <- LongSession(session)                           # Session --> long_seshdata
  #--- Mesh Data                                                    ----
  mesh <- LoadMesh(threadval$seshpath)                            # Load mesh from mesh file
  long$meshdata <- LongMesh(mesh)                                 # Mesh data --> long_meshdata
  #--- Join Data                                                    ----
  long$threaddata = long$meshdata                                 # Start with LARGEST data
  long$threaddata = LongJoin(long$threaddata, long$seshdata)      # Left join with session data
  long$threaddata <- LongAirfoil(long)                            # Special join for airfoil data
  #--- Local Mesh                                                   ----
  long <- LocalWallH(long)                                        # Average height of elements at wall
  long$threaddata <-  LongJoin(                                   # Join with threaddata
    long$threaddata, select(long$wallaveh, enum, base, aveh, ar))
  long$local <- LocalMesh(long)                                   # Determine how local the sessions elements are
  long$threaddata <- LongJoin(
    long$threaddata, select(long$local, -nnum))
  # Offset
  long$offset <- AirfoilOffset(long, totdist = 0.008, varh = TRUE)
  # ggplot() +
  #   geom_polygon(aes(x, y, group = enum), fill = NA, colour = "grey",
  #                data = long$threaddata %>% filter(!is.na(nnum), local <= 2) %>% arrange(enum, ncorner)) +
  #   geom_point(aes(x, y), alpha = 0.2, shape = 'o',
  #              data = long$threaddata %>% filter(local <= 2)) + 
  #   geom_point(aes(x, y, colour = nstep), alpha = 0.8,
  #              data = long$offset) + 
  #   # coord_fixed(xlim = c(0.55, 0.7)) +
  #   coord_fixed(xlim = c(-0.45, -0.3)) +
  #   # coord_fixed() +
  #   scale_colour_gradientn(colours = spectralpalette(6))
  #--- Airfoil Transform                                          ----
  long <- AirfoilTransform(long, localnum = 2)
  

}

# Assume:
#  Since dump files come from session files, the same
#  long$... can be used for all the dump files
#--- Dump file ----
dumpfile = LoadDump(threadval$folder, threadval$dumpfile)
long_dump <-  cbind(long_meshdata, dumpfile$flowfield)
long_localdump <- filter(long_dump, mnum %in% long_localmesh$mesh$mnum)
if (nrow(long_localdump) != nrow(long_localmesh$mesh)) warning("Not all dump nodes found") 
#--- Interpolate ----


# Function to match airfoil data to batch list -- PUT IN MAIN LOOP
BatchDumpFile <- function(dumplist, airfoillist) {
  print(dumplist$airfoil)
  airfoildata = airfoillist[dumplist$airfoil]
  # dumpfile = LoadDump(dumplist$folder, dumplist$dumpfile)   # read dump file in the main lapply loop or too much RAM
  return(list(dumplist = dumplist, airfoildata = airfoildata[[1]]))
}
# Iterate over batch list

