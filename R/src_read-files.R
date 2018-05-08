#============================
# Read Files
# Alwin Wang
#----------------------------

#--- List Session Files ----
ListSesh <- function(batchfolder) {
  # List of files in folder (should be made more robust later)
  batchlist <- list.files(batchfolder, pattern = ".sesh", recursive = TRUE)
  # Convert character to dataframe
  batchlist <- data.frame(path = unlist(strsplit(batchlist, "*.sesh"))) %>%
    mutate(seshpath = path) %>%
    separate(path, c("folder", "seshname"), sep = "/") %>%
    mutate(folder = paste0(batchfolder, "/", folder, "/"))
  # Return list of session files
  return(batchlist)
}

#--- Load Airfoil Surface Files ----
# Airfoil file should only be updated once per airfoil
# Load the wall boundary (airfoil.dat)
LoadBndry<- function(bndrypath) {
  # Session bndry name
  file = paste0(bndrypath,".dat")
  # Read bndry file
  bndry <- read.table(file)
  # Set column names
  colnames(bndry) <- c("x", "y")
  # Return bndry
  return(bndry)    # Data.frame
}
# Load the chord line
LoadChord<- function(bndrypath) {
  # Session bndry name
  file = paste0(bndrypath,"_chord.dat")
  # Read bndry file
  bndry <- read.table(file)
  # Set column names
  colnames(bndry) <- c("x", "y")
  # Return bndry
  return(bndry)    # Data.frame
}

#---- Load Mesh Files ----
# Mesh file are updated per airfoil (assuming N-order same)
# Load Mesh file generated from N order poly on mesh
LoadMesh <- function(seshpath) {
  # Session file name
  file = paste0(seshpath,".mshi")
  # Read mesh file
  mesh <- read.table(file, skip = 1)
  # Set column names
  colnames(mesh) <- c("x", "y", "enum", "jnum")
  # Return mesh
  return(mesh)    # Data.frame
}

# Load Mesh file generated from N order poly on mesh
LoadWallmsh <- function(seshpath) {
  # Session file name
  file = paste0(seshpath,".wallmsh")
  # Read mesh file
  wallmesh <- read.table(file, skip = 1)
  # Set column names
  colnames(wallmesh) <- c("x", "y")
  # Return mesh
  return(wallmesh)    # Data.frame
}

#--- Load Session File ----
# Session files are updated per session-
# Function to load keywords from the session file
LoadKeyword <- function(keyword, filelines, file) {
  # Find the line number the keyword appears in
  num <- grep(keyword[1], filelines)
  # Create a new variable with the keyword as the name
  #     could use "assign(..., envir = .GlobalEnv)" if I wanted to
  table <- read.table(
    file = file,
    skip = num[1],
    nrow = num[2]-num[1]-1,
    stringsAsFactors = FALSE)
  # Set talble colname names
  colnames(table) <- keyword[2:length(keyword)]
  # Remove any junk
  table$junk <- NULL
  return(table)     # Data.fable
}

# Function to load a list of keywords and colnames
LoadSeshFileKeywords <- function(seshpath, keywords) {
  # Session file name
  file = paste0(seshpath,".sesh")
  # Read the session file to grep lines later
  filelines <- readLines(file)
  # For each keyword load the associated data
  sessionfile <- lapply(keywords, LoadKeyword, filelines, file)
  # Set the names for the list
  names(sessionfile) <- sapply(keywords, function(key) tolower(key[1]))
  # Optional: str(sessionfile) gives short display of the contents
  # Return the session file contents
  return(sessionfile)   # List
}

# Dev temp variables
# bctext = "MOD_ALPHA_X"
# Function to read in BC equations
# Note: This cannot handle the inflow_u condition yet...
LoadSeshBCEqs <- function(seshpath, bctext, bcfuncname = NULL) {
  # Determine function name
  if (is.null(bcfuncname)) bcfuncname = bctext
  # Session file name
  file = paste0(seshpath,".sesh")
  # Read the session file to grep lines later
  filelines <- readLines(file)
  # Select line of interest
  bc <- filelines[grep(bctext, filelines)]
  # Clean up line
  bc %<>% gsub(bctext, "", .) %>%
    gsub("\t", "", .) %>%
    gsub("=", "", .) %>%
    gsub("PI", "pi", .)
  # Create the required bc function
  eval(parse(text = paste0(
    "BC_", tolower(bctext), " <- function(t) ", bc)),
    envir = .GlobalEnv)
}


#--- Load History File ----
# Load history file
LoadHist <- function(seshpath) {
  # Session file name
  file = paste0(seshpath,".his")
  # Read history file
  his <- read.table(file)
  # Set column names
  colnames(his) <- c("bnum", "t", "u", "v", "p")
  # Return his
  return(his)    # Data.frame
}

#--- Load Dump File ----
# List dump files in directory
ListDump <- function(folder, seshname) {
  # List of files in folder (should be made more robust later)
  dumplist <- list.files(folder, pattern = paste0(seshname,"-"))
  # Return file list
  return(dumplist)  # character
}

# Read single dump file
# Dev variables 
# dumpfile = dumplist[4]
LoadDump <- function(folder, dumpfile) {
  # Path to dump file
  dumppath = paste0(folder, dumpfile)
  # Read the dump file to grep lines later
  filelines <- readLines(dumppath)
  # Read time
  time <-  filelines[grep("Time", filelines)[1]] %>%
    gsub("Time", "", .) %>%
    as.numeric(.)
  # Read Kinvis
  kinvis <-  filelines[grep("Kinvis", filelines)[1]] %>%
    gsub("Kinvis", "", .) %>%
    as.numeric(.)
  # Read table
  flowfield <- read.table(
    file = dumppath,
    skip = grep("ASCII", filelines),
    stringsAsFactors = FALSE)
  colnames(flowfield) <- c("u", "v", "p", "t")
  # Return list
  return(
    list(time = time, kinvis = kinvis, flowfield = flowfield)
  )
}

