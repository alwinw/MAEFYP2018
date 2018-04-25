#============================
# Read Files
# Alwin Wang
#----------------------------

# Temp Variables
folder = "../session-files/NACA0012-AoA04/"
seshname = "RE-10000-sine-0.001-2000"
seshpath = paste0(folder, seshname)

bndry = "bndry_prf"

keywords <- list(
  c("NODES", "nnum", "x", "y", "z"),
  c("ELEMENTS", "enum", "shapetag", "n1", "n2", "n3", "n4", "junk"),
  c("SURFACES", "snum", "element", "side", "bctag", "bc", "junk"),
  c("CURVES", "cnum", "element", "side", "curvetag", "curvedata", "junk")
)
str(keywords)

#--- Load Boundary Files ----

#--- Load Session File ----
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
  return(table)     # Datatable
}

# Function to load a list of keywords and colnames
LoadSeshFile <- function(seshpath, keywords) {
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

#---- Combine Data ----

