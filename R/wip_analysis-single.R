#============================
# Analysis of Single File
# Alwin Wang
#----------------------------

#--- Structure ----
# - Intial set up
# - Run once per folder (airfoil)
# - Run once per session file (Reynolds number)
# - Run once per dump file (time step)

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
# Preprocess read files
source("src_airfoil-analysis.R")

theme_set(theme_bw())

#--- Load and Pre-process Data ----
# Location
# dir("../session-files/NACA0012-AoA04")
folder = "../session-files/NACA0012-AoA04/"
seshname = "RE-10000-sine-0.001-2000"
seshpath = paste0(folder, seshname)
bndryname = "bndry_prf"
bndrypath = paste0(folder,bndryname)
# Keywords to load
keywords <- list(
  c("NODES", "nnum", "x", "y", "z"),
  c("ELEMENTS", "enum", "shapetag", "n1", "n2", "n3", "n4", "junk"),
  c("SURFACES", "snum", "element", "side", "bctag", "bc", "junk"),
  c("CURVES", "cnum", "element", "side", "curvetag", "curvedata", "junk")
)
# Load airfoil surface
bndry <- LoadBndry(bndrypath)
chord <- LoadChord(bndrypath)

# Load mesh files
mesh <- LoadMesh(seshpath)
wallmsh <- LoadWallmsh(seshpath)
# Load session file
sessionkey <- LoadSeshFileKeywords(seshpath, keywords)
list2env(sessionkey, envir = .GlobalEnv)
LoadSeshBCEqs(seshpath, "MOD_ALPHA_X")
# Load history file
his <- LoadHist(seshpath)
# Load dump file
dumplist <- ListDump(folder, seshname)
dumpfile <- LoadDump(folder, dumplist[4])
# Pre-process Data
long_wall <- AirfoilLongWall(wallmsh)
long_wall <- AirfoilSpline(long_wall)
long_wall <- AirfoilOffset(long_wall)

# source("src_pre-processing.R")

# Dump File


# Determine acceleration field along the airfoil surface

# Determine pressure gradients

# Determine circulation/vorticity (separate branch)
