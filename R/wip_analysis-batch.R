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


#--- Aircoil Calculation ----
# Determine unique airfoil types
airfoillist <- split(batchlist, batchlist$airfoil)
airfoillist <- lapply(airfoillist, function(x) x[1,])

LoadAirfoil <- function(airfoillist) {
  # Print airfoil name
  print(airfoillist$airfoil)
  # Source required functions
  source("src_read-files.R")
  # Read the bndry file
  bndrypath = paste0(airfoillist$folder, "bndry_prf")
  print(bndrypath)
  bndry <- LoadBndry(bndrypath)
  # Read the wallmesh file
  wallmsh <- LoadWallmsh(airfoillist$seshpath)
  # Determine quantities of interest for long_wall
  long_wall <- AirfoilLongWall(wallmsh)
  long_wall <- AirfoilSpline(long_wall)
  # Return the output as a list
  return(list(airfoil = airfoillist, bndry = bndry, long_wall = long_wall))
}

airfoillist <- lapply(airfoillist, LoadAirfoil)

MatchAirfoil <- function(batchlist, airfoillist) {
  print(batchlist$airfoil)
  airfoil_data = airfoillist[batchlist$airfoil]
  return(list(batchlist = batchlist, airfoil_data = airfoil_data[[1]]))
}

threadlist <- split(batchlist, batchlist$seshpath)

temp2 <- lapply(threadlist, MatchAirfoil, airfoillist)

str(temp2[[1]])

# Start the cluster
cl <- makeCluster(detectCores())
# Export objects into the cluster
clusterExport(cl, c("airfoillist"))
# Function to be applied across the list



temp <- pbapply(airfoillist, 1, tempfunc, cl = cl)


# Stop cluster
stopCluster(cl)
