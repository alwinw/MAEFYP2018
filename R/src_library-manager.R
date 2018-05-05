#============================
# Library Manager
# Alwin Wang
#----------------------------

#--- Recommendations ----
# R version: 64-bit https://cran.r-project.org/bin/windows/base/
# RStudio: Free Versions https://www.rstudio.com/products/rstudio/download/
# Use RStudio projects to manage code and workspace

#--- List of  Required Packages ----
# List of packages required
packages <- c(
  "stats",      # For statistical analysis, Part of R Core
  "pracma",     # https://www.rdocumentation.org/packages/pracma/versions/1.9.9
  "ggplot2",    # Great for visualising data, by Adj. Prof Hadley Wickham
  "dplyr",      # Great for manipulating data, by Adj. Prof Hadley Wickham
  "tidyr",      # Great for tidying data, by Adj. Prof Hadley Wickham
  "stringr",    # Great for woking with strings, by Adj. Prof Hadley Wickham
  "lubridate",  # Great for working with dates, by Adj. Prof Hadley Wickham
  "readxl",     # Great for reading .xls, .xlsx, by Adj. Prof Hadley Wickham 
  "devtools",   # Great for developing, by Adj. Prof Hadley Wickham 
  "rstudioapi", # For setting working directory
  "RColorBrewer"  # For ggplot2 colours 
)

#--- Install Missing Packages ----
# install.packages("tidyverse")
# Check if any packages need to be installed
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
# Install missing packages
if(length(new.packages)) install.packages(new.packages)

#--- Load Libraries ----
invisible(lapply(packages, library, character.only = TRUE))
# Clean up
rm(packages, new.packages)