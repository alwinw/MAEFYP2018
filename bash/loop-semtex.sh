#!/bin/bash

#============================
# Loop over session files and solve using semtex dns
# Alwin Wang
#----------------------------

#--- Inputs ----
# usage="$(loop-semtex "$0") [-h] -- script to loop over session files in a folder"
# write getops case handing to display for -h
# Input 1: folder location

#--- Loop over files in folder ----
# bash                                # Enter bash terminal if not already
# echo "In bash"

#--- Working Directories ----
# Save old working dir
oldpwd=$(pwd)
echo $oldpwd
# Navigate to Folder location
folder=../session-files/
cd $folder

#--- Loop Over Airfoils in Folder ----
# Airfoil folders
airfoils=*/
# Loop over airfoils
for airfoil in airfoils; do
  echo
  echo "#============================"
  # Echo airfoil name
  echo " Found $airfoil"
  # Navigate into folder
  cd $airfoils
  echo " > Entering $(pwd)"
  #--- Loop Over Files in Folder ----
  # Session files
  sessionfiles=*.sesh
  # Loop over session files
  for sessionfile in $sessionfiles; do
    echo "#----------------------------"    
    # Remove extension
    f="${sessionfile%.*}"                     && echo "   Found $f"
  done
  # Return to previous folder
  cd ..
done

#--- Return to previous wd ----
cd $oldpwd
