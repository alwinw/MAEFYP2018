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
source("src_plot-output.R")                                     # Plots to output and save
# Custom ggplot2 setup
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
# Output Location
saveplot = "../output-plot"
# savedata = "Output_Data"
if (!dir.exists(saveplot)) dir.create(saveplot, recursive = TRUE)
# if (!file.info(savedata)$isdir) dir.create(savedata, recursive = TRUE)
# logfile = paste0(format(Sys.time(), "%Y-%m-%dT%H.%M.%S"), ".txt")

#--- List of Session FIles                                        ----
batchfolder = "../session-files"                                # Path to session files
batchlist <- ListSesh(batchfolder)                              # List sessionfiles in folder

# Check if all file types exist, if not then call bash script

#--- Airfoil Calculation                                          ----
# Function to load airfoil data
BatchLoadAirfoil <- function(airfoilval) {                      # airfoilval = airfoillist[[1]]
  source("src_library-manager.R")                                 # Call libraries and install missing ones
  source("src_numerical-methods.R")                               # Load custom numerical methods
  source("src_load-files.R")                                      # Load data
  source("src_airfoil-analysis.R")                                # Airfoil files
  source("src_plot-output.R")                                     # Plots to output and save
  #--- Boundary Data                                                ----
  bndrypath = paste0(airfoilval$folder, "bndry_prf")              # Path to bndry file
  bndry <- LoadBndry(bndrypath)                                   # Read the bndry file
  #--- Wall Mesh Data                                               ----
  wallmsh <- LoadWallmsh(airfoilval$seshpath)                     # Read the wallmesh file
  long_wall <- AirfoilLongWall(wallmsh)                           # Airfoil data --> long_wall
  long_wall <- AirfoilSpline(long_wall)                           # Determine spline distance
  # Return the output as a list
  return(list(
    airfoil = airfoilval, 
    bndry = bndry, 
    long_wall = long_wall))
    # , unixy_wallmsh = unixy_wallmsh))
}

# Iterate over list of airfoils list to load data
airfoillist <- split(batchlist, batchlist$airfoil)              # Determine unique airfoil types
airfoillist <- lapply(airfoillist, function(x) x[1,])           # Take only the first row in each
# airfoillist <- lapply(airfoillist, BatchLoadAirfoil)
cl <- makeCluster(detectCores())                                # Start the cluster
clusterExport(cl, c("airfoillist", "BatchLoadAirfoil"))         # Export objects into the cluster
airfoillist <- pblapply(airfoillist, BatchLoadAirfoil, cl = cl) # Load airfoil data from each airfoil
stopCluster(cl)                                                 # Stop the cluster

#--- Session and Mesh Calculation                                 ----
# Function to load mesh data                                    # meshval = meshlist[[1]]
BatchLoadMesh <- function(meshval, airfoillist) {
  source("src_library-manager.R")                                 # Call libraries and install missing ones
  source("src_numerical-methods.R")                               # Load custom numerical methods
  source("src_load-files.R")                                      # Load data
  source("src_airfoil-analysis.R")                                # Airfoil files
  source("src_vorticity-generation.R")                            # Vorticity Generation
  source("src_plot-output.R")                                     # Plots to output and save
  long <- list()                                                  # Create a list of long format data
  #--- Airfoil Data                                                 ----
  airfoildata = airfoillist[[meshval$airfoil]]                    # Collect airfoil data
  long$walldata = airfoildata$long_wall                           # Airfoil --> long_walldata
  #--- Session Data                                                 ----
  keywords <- list(                                               # Keywords in session to read
    c("NODES", "nnum", "x", "y", "z"),
    c("ELEMENTS", "enum", "shapetag", "n1", "n2", "n3", "n4", "junk"),
    c("SURFACES", "snum", "element", "side", "bctag", "bc", "junk"),
    c("CURVES", "cnum", "element", "side", "curvetag", "curvedata", "junk"))
  session <- LoadSeshFileKeywords(meshval$seshpath, keywords)     # Read keywords from session file
  long$seshdata <- LongSession(session)                           # Session --> long_seshdata
  #--- Mesh Data                                                    ----
  mesh <- LoadMesh(meshval$seshpath)                              # Load mesh from mesh file
  long$meshdata <- LongMesh(mesh)                                 # Mesh data --> long_meshdata
  #--- Join Data                                                    ----
  long$threaddata = long$meshdata                                 # Start with LARGEST data
  long$threaddata = LongJoin(long$threaddata, long$seshdata)      # Left join with session data
  long$threaddata <- LongAirfoil(long)                            # Special join for airfoil data
  long$threaddata$wall = !is.na(long$threaddata$wnum)
  long$threaddata$seshnode = !is.na(long$threaddata$nnum)
  #--- Local Mesh                                                   ----
  long <- LocalWallH(long)                                        # Average height of elements at wall
  long$threaddata <-  LongJoin(                                   # Join with threaddata
    long$threaddata, select(long$wallaveh, enum, base, aveh, ar))
  long$local <- LocalMesh(long)                                   # Determine how local the sessions elements are
  long$threaddata <- LongJoin(
    long$threaddata, select(long$local, -nnum))
  # Offset
  long$offset <- AirfoilOffset(long,                              # Create df of offset points from surface
                               totdist = 0.008, varh = TRUE)
  long <- AirfoilOffsetEnum(long)                                 # Update enum values of the offset points
  #--- Airfoil Transform                                          ----
  long <- AirfoilTransform(long, localnum = 2)
  # Do some plots and save them to png
  #--- Clean up and Return                                        ----
  # Clean up
  output <- long[c("walldata", "threaddata", "offset")]
  return(output)
}

# Load session info e.g. tokenwords = list("N_P", "N_Z")
tokenwords = list("N_P")                                        # Find unique bndry and knot N values
batchlist <- batchlist %>%
  rowwise() %>%
  do(data.frame(., LoadSeshTokenWords(.$seshpath, tokenwords),
                stringsAsFactors = FALSE)) %>%
  mutate(ID = paste(airfoil, 
                    paste0(tokenword, tokenvalue), sep = "-"))  # This CANNOT handle multiple token words...
meshlist <- split(batchlist, batchlist$ID)                      # Unique airfoil and N_P
meshlist <- lapply(meshlist, function(x) x[1,])                 # Take only the first row in each
# meshlist <- lapply(meshlist, BatchLoadMesh)
cl <- makeCluster(detectCores())                                # Start the cluster
clusterExport(cl, c("meshlist", "BatchLoadMesh"))               # Export objects into the cluster
meshlist <- pblapply(meshlist, BatchLoadMesh,                   # Load airfoil data from each airfoil
                     airfoillist, cl = cl)
stopCluster(cl)                                                 # Stop the cluster

#--- Dump File List                                               ----
# Function to load and process dump files
BatchLoadDump <- function(dumpval, meshlistin) {              # dumpval = dumplist[[5]]; meshlistin = meshlist
  source("src_library-manager.R")                                 # Call libraries and install missing ones
  source("src_numerical-methods.R")                               # Load custom numerical methods
  source("src_load-files.R")                                      # Load data
  source("src_airfoil-analysis.R")                                # Airfoil files
  source("src_vorticity-generation.R")                            # Vorticity Generation
  source("src_plot-output.R")                                     # Plots to output and save
  #---  Long Mesh Data                                              ----
  long <- meshlistin[[dumpval$ID]]                                # Collect airfoil data
  rm(meshlistin)                                                  # Reduce memory required
  #--- Dump Data                                                    ----
  dump <- LoadDump(dumpval$folder, dumpval$dumpfile)              # Load dump file as list
  dump$threaddata <- LongDump(dump, long)
  #--- Acceleration                                                 ----
  LoadSeshBCEqs(dumpval$seshpath, "MOD_ALPHA_X")                  # Load BC Equation from session file
  dump$acceleration = BC_mod_alpha_x(dump$time)                   # Instantaenous acceleration
  dump$accel = DumpAcceleration(dump, long)                       # Calculate acceleration components
  dump$threaddata = LongJoin(dump$threaddata, dump$accel)         # Join with rest of data
  #--- Pressure                                                     ----
  dump$pres = DumpPressureStream(dump, long)                      # Calculate pressure gradient dp/ds
  dump$threaddata = LongJoin(dump$threaddata, dump$pres)          # Join with rest of data
  #--- Vorticity Interpolation                                      ----
  dump$offset = long$offset                                       # Initialise the offset for interpolation
  # dump <- DumpVortTransformed(dump, localval = 2, var = "t")      # Interpolate based on stream, norm cs
  dump <- DumpVortElements(dump, localval = 2, var = "t")         # Interpolate using each element
  # Derivatives
  # dump$offset <- FiniteDiff(dump$offset, "t_trans")               # Calculate derivative using stream, norm cs
  dump$offset <- FiniteDiff(dump$offset, "t_enum")                # Calculate derivative using each element
  #--- Plot Output                                                  ----
  # Vertical Lines
  plot_vline <- data.frame(                                       # LE, TE, LE limits for vertical lines
    te_up = min(long$walldata$s),
    le = long$walldata[long$walldata$theta == pi,]$s[1],
    te_lo = max(long$walldata$s))
  plot_surf <- data.frame(                                        # Labels for surfaces
    x = c(plot_vline$te_up,
          (plot_vline$te_up + plot_vline$le)/2,
          plot_vline$le,
          (plot_vline$le + plot_vline$te_lo)/2,
          plot_vline$te_lo),
    y = rep(-20, 5),
    labels = c("TE", "Upper Surface", "LE", "Lower Surface", "TE"))
  plot_xbreaks <-                                                 # Breaks in x axis
    c(seq(plot_vline$te_up, plot_vline$le, length.out = 3)[1:2],
      seq(plot_vline$le, plot_vline$te_lo, length.out = 3))
  # Title
  plot_title <- paste(
    dumpval$airfoil,
    paste("Kinematic Viscosity:", format(dump$kinvis, scientific = FALSE)),
    paste("Time:", sprintf("%4.2f", dump$time)),
    paste("Acceleration:", sprintf("%+6.3f", dump$acceleration)),
    sep = "\n")
  plot_filename <- paste(
    dumpval$ID,
    paste0("v", format(dump$kinvis, scientific = TRUE)),
    paste0("t", sprintf("%05.3f", dump$time)),
    paste0("a", sprintf("%+07.3f", dump$acceleration)),
    sep = "_")
  # Plot
  plot_NS <- ggplot(dump$threaddata %>%                           # Initiate plot with threaddata
                   filter(wall) %>% arrange(s), aes(s)) +
    geom_vline(xintercept = as.numeric(plot_vline),               # Vertical lines for LE, TE, LE
               colour = "grey", linetype = "dashed") +
    geom_label(aes(x, y, label = labels), plot_surf,              # Surface labels
              colour = "grey") +
    geom_path(aes(y = -accels, colour = "dU/dt")) +               # Acceleration terms
    geom_path(aes(y = dpds, colour = "dp/ds")) +                  # Pressure field
    geom_path(aes(y = - accels + dpds, colour = "LHS"),           # LHS, acceleration + pressure
              linetype = "dashed") +
    geom_point(aes(s, t_enum_diff*dump$kinvis, colour = "RHS"),   # RHS, v * dw/dz
               dump$offset, shape = "O") +
    xlab("s") + 
    scale_x_continuous(breaks = plot_xbreaks, 
                       labels = function(x) sprintf("%.2f", x)) +
    ylab(NULL) + ylim(c(-20, 40)) +
    scale_color_manual(
      name = "Legend",
      values = c("dp/ds" = "red", "dU/dt" = "blue", "LHS" = "purple", "RHS" = "purple"),
      labels = c(
        expression(frac(1, rho)~frac(partialdiff*p, partialdiff*s)), 
        expression(frac(partialdiff*U, partialdiff*t)), 
        expression(
          bgroup("(",
            frac(1, rho)~frac(partialdiff*p, partialdiff*s) +
              frac(partialdiff*U, partialdiff*t), ")")),
        expression(nu~frac(partialdiff*omega, partialdiff*z))),
      guide = guide_legend(
        override.aes = list(
          linetype = c(rep("solid", 2), "dashed", "blank"),
          shape = c(rep(NA, 3), "O")))) +
    theme(legend.key.size = unit(2.25, "lines"),
          legend.text.align = 0.5,
          legend.direction = "vertical", 
          legend.position = "right",
          legend.background = element_rect(colour = "black", size = 0.3)) +
    ggtitle(plot_title)
  
  # Should save to indivial folders based on ID
  ggsave(paste0(plot_filename, ".png"), plot_NS, path = "../output-plot",
         width = 8, height = 6)
  # Return the output
  return(dump)
}
# Process batch list
batchlist <- batchlist %>%                                      # Determine dump list
  rowwise() %>% 
  do(data.frame(., dumpfile = ListDump(.$folder, .$seshname),
                stringsAsFactors = FALSE)) %>%
  mutate(dumppath = paste0(folder, dumpfile))
dumplist <- split(batchlist, batchlist$dumppath)                # Create thread list
# Testing, shortern dump list
# dumplist <- dumplist[1:8]
# dumplist <- lapply(dumplist, BatchLoadDump)
cl <- makeCluster(detectCores())                                # Start the cluster
clusterExport(cl, c("dumplist", "BatchLoadDump"))               # Export objects into the cluster
dumplist <- pblapply(dumplist, BatchLoadDump,                   # Load airfoil data from each airfoil
                     meshlist, cl = cl)
stopCluster(cl)   
