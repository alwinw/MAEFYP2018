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
source(paste0(srcpath, "src_helper-functions.R"))               # Smaller functions used
# Additional scripts here
# ggplot2 setup (consider moving to a separate script)
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 

#--- * Required Paths                                             ----
if (FALSE) {
  #--- > Single File Analysis                                       ----
  saveplot   = "images/"
  airfoil    = "NACA0012r"
  folderpath = "results/"
  seshpath   = "NACA0012r_kutta"
  dumppath   = "NACA0012r_kutta-249.dump"
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
  outp_dump           <- list(BatchLoadDump(data_dump, outp_mesh, 
                                            plot = c("images/", "airfoil", "TEstream"), 
                                            outp = "te", srcpath))
  names(outp_dump)    <- data_dump$dumppath
  
} else {
  # Batch Test Case
  # Batch list to process
  batchfolder = "results"
  df_batch <- ListSesh(batchfolder, c("airfoil"))
  # List of unique airfoils & N_P
  li_airfoil   <- split (df_batch, df_batch$airfoil)
  li_airfoil   <- lapply(li_airfoil, function(df) df[1,])
  outp_airfoil <- pblapply(li_airfoil, BatchLoadAirfoil, srcpath)
  # List of meshs (should turn this into a function!)
  tokenwords = list("N_P")                                        # Find unique bndry and knot N values
  df_mesh <- df_batch %>%
    rowwise() %>%
    do(data.frame(
      ., LoadSeshTokenWords(.$seshpath, tokenwords), stringsAsFactors = FALSE)) %>%
    filter(!is.na(tokenvalue)) %>% 
    mutate(                                                       # May need to adjust later
      ID = paste0(airfoil))      
  li_mesh   <- split (df_mesh, df_mesh$ID)                        # Unique airfoil and N_P
  li_mesh   <- lapply(li_mesh, function(x) x[1,])                 # Take only the first row in each
  outp_mesh <- pblapply(li_mesh, BatchLoadMesh, outp_airfoil, srcpath)
  # List of dumps
  df_dump <- df_mesh %>% 
    rowwise() %>% 
    do(data.frame(., dumpfile = ListDump(.$folder, .$seshname),
                  stringsAsFactors = FALSE)) %>%
    mutate(dumppath = paste0(folder, dumpfile))
  li_dump   <- split(df_dump, df_dump$dumppath)
  # Use cluster
  cl <- makeCluster(detectCores() - 1)
  # clusterExport(cl, c("outp_mesh", "srcpath"))
  outp_dump <- pblapply(li_dump, BatchLoadDump, outp_mesh, 
                        plot = c("images/", "airfoil", "TEstream"), 
                        outp = "te", srcpath,
                        cl = cl)
  stopCluster(cl)
}
