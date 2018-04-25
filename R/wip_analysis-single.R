#============================
# Analysis of Single File
# Alwin Wang
#----------------------------

#--- Set Working Directory ----
# Use rstudioapi to get saved location of this file
# and set it as the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Requres devtools, rstudioapi

#--- Source Required Scripts ----
# Call libraries and install missing ones
source("src_library-manager.R")




#--- Load File to Analyse ----
# dir("../session-files/NACA0012-AoA04")
folder = "../session-files/NACA0012-AoA04/"
file = "RE-10000-sine-0.001-2000-02.dump"

temp = readLines(paste0(folder,file))
