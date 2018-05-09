#============================
# Batch Analysis
# Alwin Wang
#----------------------------

#--- Set Working Directory ----
# Use rstudioapi to get saved location of this file
# and set it as the working directory
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
BatchLoadAirfoil <- function(airfoillist) {
  # Source required functions
  source("src_read-files.R")
  # Read the bndry file
  bndrypath = paste0(airfoillist$folder, "bndry_prf")
  bndry <- LoadBndry(bndrypath)
  # Read the wallmesh file
  wallmsh <- LoadWallmsh(airfoillist$seshpath)
  # Determine quantities of interest for long_wall
  long_wall <- AirfoilLongWall(wallmsh)
  long_wall <- AirfoilSpline(long_wall)
  # Return the output as a list
  return(list(airfoil = airfoillist, bndry = bndry, long_wall = long_wall))
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
  # Collect airfoil data
  airfoildata = airfoillist[[threadval$airfoil]]
  # Collect dump file
  dumpfile = LoadDump(threadval$folder, threadval$dumpfile)
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
