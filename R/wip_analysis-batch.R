#============================
# Batch Analysis
# Alwin Wang
#----------------------------

#--- Set Working Directory ----
# Use rstudioapi to get saved location of this file
# https://rviews.rstudio.com/2016/11/11/easy-tricks-you-mightve-missed/
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Requres devtools, rstudioapi

#--- Source Required Scripts ----
# Call libraries and install missing ones
source("src_library-manager.R")
# Load custom numerical methods
source("src_numerical-methods.R")
# Read data
source("src_read-files.R")
# Airfoil files
source("src_airfoil-analysis.R")
# Vorticity Generation
source("src_vorticity-generation.R")
# Custom setup
theme_set(theme_bw())

#--- Load List of Session FIles ----
# Compile a list of the various session files
batchfolder = "../session-files"
batchlist <- ListSesh(batchfolder, listout = FALSE)

# Check if all file types exist, if not then call bash script

#--- Output Location ----
# saveplot = "Output_Plot"
# savedata = "Output_Data"
# if (!file.info(saveplot)$isdir) dir.create(saveplot, recursive = TRUE)
# if (!file.info(savedata)$isdir) dir.create(savedata, recursive = TRUE)
# logfile = paste0(format(Sys.time(), "%Y-%m-%dT%H.%M.%S"), ".txt")


#--- Airfoil Calculation ----
# Determine unique airfoil types
airfoillist <- split(batchlist, batchlist$airfoil)
airfoillist <- lapply(airfoillist, function(x) x[1,])
# Function to load airfoil data
# airfoilval = airfoillist[[1]]
BatchLoadAirfoil <- function(airfoilval) {
  # Source required functions
  source("src_read-files.R")
  # Read the bndry file
  bndrypath = paste0(airfoilval$folder, "bndry_prf")
  bndry <- LoadBndry(bndrypath)
  # Read the wallmesh file
  wallmsh <- LoadWallmsh(airfoilval$seshpath)
  # Determine quantities of interest for long_wall
  long_wall <- AirfoilLongWall(wallmsh)
  long_wall <- AirfoilSpline(long_wall)
  # To join over (x, y), these duplicate coordinates need to be removed!
  unixy_wallmsh <- long_wall[!duplicated(select(long_wall, x, y)),]
  # Assume the mesh can change depending on n x n, so do not read here
  # Return the output as a list
  return(list(airfoil = airfoilval, bndry = bndry, long_wall = long_wall, unixy_wallmsh = unixy_wallmsh))
}
# Iterate over list of airfoils list to load data
airfoillist <- lapply(airfoillist, BatchLoadAirfoil)

# Start the cluster
cl <- makeCluster(detectCores())
# Export objects into the cluster
clusterExport(cl, c("airfoillist"))
# Function to be applied across the list
#    pblappy here
# Stop cluster
stopCluster(cl)


#---- Dump File List ----
# Determine dump list
dumplist <- batchlist %>%
  rowwise() %>% 
  do(data.frame(., dumpfile = ListDump(.$folder, .$seshname))) %>%
  mutate(dumppath = paste0(folder, "/", dumpfile))
# Create thread list
threadlist <- split(dumplist, dumplist$dumppath)  # List

#--- Main Apply Function ----
# Function for each each threadval = threadlist[[1]]
BatchThread <- function(threadval, airfoillist) {
  print(threadval$dumppath)
  #---- Airfoil Data ----
  # Collect airfoil data
  airfoildata = airfoillist[[threadval$airfoil]]
  #   airfoil data: airfoil, bndry, long_wall
  long_wall = airfoildata$long_wall
  #--- Session Data ----
  # Keywords in session to read
  keywords <- list(
    c("NODES", "nnum", "x", "y", "z"),
    c("ELEMENTS", "enum", "shapetag", "n1", "n2", "n3", "n4", "junk"),
    c("SURFACES", "snum", "element", "side", "bctag", "bc", "junk"),
    c("CURVES", "cnum", "element", "side", "curvetag", "curvedata", "junk"))
  # Consider combining all the following items into lists e.g. session <- list(), mesh <- list()
  # Read keywords from session file
  session <- LoadSeshFileKeywords(threadval$seshpath, keywords)
  # Load BC Equation from session file
  LoadSeshBCEqs(threadval$seshpath, "MOD_ALPHA_X")
  # Session mesh data in long format
  long_seshdata <- SessionMesh(session)
  #---Mesh Data ----
  # Load mesh
  mesh <- LoadMesh(threadval$seshpath)
  # mesh data into long format
  long_meshdata <- MeshLong(mesh, airfoildata)
  #--- Local Mesh ----
  # local mesh around the airfoil
  long_localmesh <- MeshLocal(long_seshdata, long_meshdata) 
  # Airfoil offset
  airfoildata$offset <-  AirfoilOffset(long_wall, totdist = long_localmesh$localave$mean, nsteps = 5) 
  #--- Dump file ----
  dumpfile = LoadDump(threadval$folder, threadval$dumpfile)
  long_dump <-  cbind(long_meshdata, dumpfile$flowfield)
  long_localdump <- filter(long_dump, mnum %in% long_localmesh$mesh$mnum)
  if (nrow(long_localdump) != nrow(long_localmesh$mesh)) warning("Not all dump nodes found") 
  #--- Interpolate ----
}



# Function to match airfoil data to batch list -- PUT IN MAIN LOOP
BatchDumpFile <- function(dumplist, airfoillist) {
  print(dumplist$airfoil)
  airfoildata = airfoillist[dumplist$airfoil]
  # dumpfile = LoadDump(dumplist$folder, dumplist$dumpfile)   # read dump file in the main lapply loop or too much RAM
  return(list(dumplist = dumplist, airfoildata = airfoildata[[1]]))
}
# Iterate over batch list
threadlist <- split(dumplist, dumplist$seshpath)
threadlist <- lapply(threadlist, MatchAirfoil, airfoillist)

str(threadlist[[1]])
