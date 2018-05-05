#============================
# Analysis of Single File
# Alwin Wang
#----------------------------

theme_set(theme_bw())

#--- Set Working Directory ----
# Use rstudioapi to get saved location of this file
# and set it as the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Requres devtools, rstudioapi

#--- Source Required Scripts ----
# Call libraries and install missing ones
source("src_library-manager.R")
# Read data
source("src_read-files.R")
# Preprocess read files

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
session <- LoadSeshFile(seshpath, keywords)
list2env(session, envir = .GlobalEnv)

