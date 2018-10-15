#============================#
# Analysis of Kutta Condition
# Alwin Wang
#----------------------------#

#--- Set Up                                                       ----
# Use rstudioapi to get saved location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
#--- * Scripts                                                    ----
# Source Required Scripts
srcpath = "../../R/"
source(paste0(srcpath, "src_library-manager.R"))                # Call libraries and install missing ones
source(paste0(srcpath, "src_batch-functions.R"))                # Batch functions used
# source(paste0(srcpath, "src_helper-functions.R"))               # Smaller functions used
# Additional scripts here
# ggplot2 setup (consider moving to a separate script)
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 

#--- * Required Paths                                             ----
if (TRUE) {
  #--- > Single File Analysis                                       ----
  saveplot   = "images/"
  airfoil    = "NACA0012r"
  folderpath = "results/"
  seshpath   = "NACA0012r_kutta"
  dumppath   = "NACA0012r_kutta-180.dump"
  #--- > Required input dataframes                                  ----
  data_airfoil <- data.frame(
    airfoil  = airfoil,
    seshname = seshpath,
    folder   = folderpath,
    seshpath = paste0(folderpath, seshpath),
    stringsAsFactors = FALSE)
  data_mesh <- data.frame(
    data_airfoil,
    tokenword  = "N_P",
    tokenvalue = 8,
    ID         = "NACA0012r-N_P8",
    stringsAsFactors = FALSE)
  data_dump <- data.frame(
    data_mesh,
    dumpfile = dumppath,
    dumppath = paste0(folderpath, dumppath),
    stringsAsFactors = FALSE)
  rm(airfoil, folderpath, seshpath, dumppath)
  #--- > Call Batch Functions                                       ----
  plot = 0
  # Airfoil
  outp_airfoil        <- list(BatchLoadAirfoil(data_airfoil, srcpath))
  names(outp_airfoil) <- data_airfoil$airfoil
  # Mesh
  outp_mesh           <- list(BatchLoadMesh(data_mesh, outp_airfoil, srcpath))
  names(outp_mesh)    <- data_mesh$ID
  # Dump
  outp_dump           <- list(BatchLoadDump(data_dump, outp_mesh, plot, srcpath))
  names(outp_dump)    <- data_dump$dumppath
  
} else {
  # Batch Test Case
}
