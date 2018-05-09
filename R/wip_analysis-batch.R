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
# Custom setup
theme_set(theme_bw())

#--- Load List of Session FIles ----
# Compile a list of the various session files
batchfolder = "../session-files"
batchlist <- ListSesh(batchfolder, listout = TRUE)

# Check if all file types exist, if not then call bash script

#--- Output Location ----
# saveplot = "Output_Plot"
# savedata = "Output_Data"
# if (!file.info(saveplot)$isdir) dir.create(saveplot, recursive = TRUE)
# if (!file.info(savedata)$isdir) dir.create(savedata, recursive = TRUE)
# logfile = paste0(format(Sys.time(), "%Y-%m-%dT%H.%M.%S"), ".txt")


#--- Aircoil Calculation ----
# Determine unique airfoil types
airfoillist <- lapply(batchlist, function(x) x[1,])

tempfunc2 <- function(airfoillist) {
  # Print airfoil name
  print(airfoillist$airfoil)
  # Source required functions
  source("src_read-files.R")
  bndrypath = paste0(airfoillist$folder, "bndry_prf")
  print(bndrypath)
  bndry <- LoadBndry(bndrypath)
  return(bndry)
}

temp <- lapply(airfoillist, tempfunc2)



# Start the cluster
cl <- makeCluster(detectCores())
# Export objects into the cluster
clusterExport(cl, c("airfoillist"))
# Function to be applied across the list

tempfunc <- function(airfoillist) {
  # Since this is in an apply function, atomic not recursive object. So, use "getElement" or [ ]
  source("src_read-files.R")
  bndrypath = paste0(airfoillist["folder"], "bndry_prf")
  print(bndrypath)
  bndry <- LoadBndry(bndrypath)
  return(bndry)
}

temp <- pbapply(airfoillist, 1, tempfunc, cl = cl)


# Stop cluster
stopCluster(cl)
