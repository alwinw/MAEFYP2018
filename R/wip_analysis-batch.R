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
  #--- Session Data ----
  # Keywords in session to read
  keywords <- list(
    c("NODES", "nnum", "x", "y", "z"),
    c("ELEMENTS", "enum", "shapetag", "n1", "n2", "n3", "n4", "junk"),
    c("SURFACES", "snum", "element", "side", "bctag", "bc", "junk"),
    c("CURVES", "cnum", "element", "side", "curvetag", "curvedata", "junk"))
  # Read keywords from session file
  session <- LoadSeshFileKeywords(threadval$seshpath, keywords)
  # Load BC Equation from session file
  LoadSeshBCEqs(threadval$seshpath, "MOD_ALPHA_X")
  # Elements
  # Check that all elements are quadrangles
  if (as.logical(sum(session$elements$shapetag != "<Q>"))) warning("Not all elements are Quadrangles") 
  # Manipulate data
  long_seshdata <- session$elements %>%
    # Remove shape tag since checked shape already
    select(-shapetag) %>%
    # Gather nodes 1-5 into a single column
    gather(ncorner, nnum, -enum) %>%
    # Combine with nodes to get (x,y) for each element corner
    left_join(., session$nodes, "nnum") %>%
    # Order data
    group_by(enum) %>%
    # Add columns for labels
    mutate(
      elabx = mean(x),
      elaby = mean(y)) %>%
    # Close the loop (without effecting mean)
    rbind(., mutate(filter(., ncorner=="n1"), ncorner = "n5")) %>%
    arrange(enum, ncorner) %>% 
    # Calculate area
    mutate(area = x*lead(y) - lead(x)*y) %>%
    mutate(area = 1/2*abs(sum(ifelse(is.na(area), 0, area)))) %>%
    # Join with wall mesh
    left_join(., airfoildata$unixy_wallmsh, by = c("x", "y")) %>%
    mutate(wall = !is.na(wnum)) %>%
    # Remove extra ncorner
    filter(ncorner != "n5")
  #---Mesh Data ----
  # Load mesh
  mesh <- LoadMesh(threadval$seshpath)
  # Determine which mesh data belong to which
  mesh$mnum = 1:nrow(mesh)
  # Join with unique (x, y) for temp_wallmsh
  long_meshdata <- left_join(mesh, airfoildata$unixy_wallmsh, by = c("x", "y"))
  long_meshdata$wall = !is.na(long_meshdata$wnum)
  # Check the number of rows has not changed
  if (nrow(long_meshdata) != nrow(mesh)) {
    warning("Number of nodes has changed!")}
  # Check that all wall mesh nodes were found
  if (nrow(airfoildata$unixy_wallmsh) != 
      nrow(long_meshdata %>% filter(wall) %>% select(wnum) %>% unique())) {
    warning("Not all wall mesh nodes found")}
  # Determine which mesh nodes are session node numbers
  long_meshdata <- left_join(long_meshdata, long_seshdata, by = c("enum", "x", "y"))
  long_meshdata$node = !is.na(long_meshdata$nnum)
  if (nrow(filter(long_seshdata)) !=
      nrow(long_meshdata %>% select(enum, ncorner) %>% filter(!is.na(ncorner)) %>% unique())) {
    warning("Not all mesh nodes found")}
  # Determine mesh offset
  pres_h <- long_seshdata %>% 
    # Remove non session node points, after original elements only
    group_by(enum) %>%
    # Count the number of nodes per element on the surface
    add_count(wall) %>%
    filter(n == 2, wall) %>%
    # Find the average 'height' of the rectangle with base on surface
    mutate(base = EucDist(x, y),
           aveh = area/base,
           ar = base/aveh) %>%
    # Filter results to one per element
    filter(!is.na(base))
  median(pres_h$aveh)
  min(pres_h$aveh)
  
  # Inital nodes and elements
  pres_mesh <- filter(long_seshdata, wall)
  pres_mesh <- list(
    nodes = unique(pres_mesh$nnum),
    elements = unique(pres_mesh$enum))
  # Plot
  pres_mesh$mesh <- filter(long_meshdata, enum %in% pres_mesh$elements)
  ggplot() + 
    geom_point(aes(x, y), pres_mesh$mesh, alpha = 0.2) + 
    geom_polygon(aes(x, y, group = enum), 
                 filter(pres_mesh$mesh, node) %>% arrange(ncorner), fill = NA, colour = "black") + 
    coord_fixed()
  
  # Next ndoes and elements
  pres_mesh$nodes =unique(
    long_seshdata[long_seshdata$enum %in% pres_mesh$elements,]$nnum)
  pres_mesh$elements =unique(
    long_seshdata[long_seshdata$nnum %in% pres_mesh$nodes,]$enum)
  # Plot
  pres_mesh$mesh <- filter(long_meshdata, enum %in% pres_mesh$elements)
  ggplot() + 
    geom_point(aes(x, y), pres_mesh$mesh, alpha = 0.2) + 
    geom_polygon(aes(x, y, group = enum), 
                 filter(pres_mesh$mesh, node) %>% arrange(ncorner), fill = NA, colour = "black") + 
    coord_fixed()
  
  #--- Dump file ----
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
