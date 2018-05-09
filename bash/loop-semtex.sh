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
for airfoil in $airfoils; do
  echo
  echo "#============================"
  # Echo airfoil name
  echo " Found $airfoil"
  # Navigate into folder
  cd $airfoil
  echo " > Entering $(pwd)"
  #--- Loop Over Files in Folder ----
  # Session files
  sessionfiles=*.sesh
  # Loop over session files
  for sessionfile in $sessionfiles; do
    echo "#----------------------------"    
    # Remove extension
    f="${sessionfile%.*}"                               && echo "   Found $f"             &&
      # Copy to extensionless file
      cp $sessionfile "$f"                              && echo "   > Copy Made"          &&
      # Generate the mesh
      meshpr $f > "$f.msh"                              && echo "   > Meshpr Finished"    &&
      # Add extra mesh info
      meshpr -i $f > "$f.mshi"                          && echo "   > Meshpr -i Finished" &&
      # Enumerate
      enumerate $f > "$f.num"                           && echo "   > Enumerate Finished" &&
      # Compare
      compare $f > "$f.rst"                             && echo "   > Compare Finished"   &&
      # Generate wall mesh
      wallmesh $f "$f.msh" > "$f.wallmsh"               && echo "   > Wall Mesh Finished" &&
      # dns (verbosity -v turned off);
      # use > /dev/null to hide stdout and 2>&1 to hide stderror as well
      dns $f | grep "Divergence Energy:"                && echo "   > DNS Finished"       &&
      # Add vorticity field t
      # Note: new field file CANNOT be the same as the original one (else 0)
      addfield -v -s $f "$f.fld" > "$f.vfld"            && echo "   > Vorticity Added"    &&
      # Convert to ascii
      convert "$f.vfld" > "$f.flddump"                  && echo "   > Convert Finished"   &&
      # Split the output file
      csplit -z "$f.flddump" /Session/ {*} >/dev/null   && echo "   > Split Finished"     &&
      # Remove unneded dump file
      rm "$f.flddump"                                   && echo "   > ASCII dump removed" &&
      # Rename split files
      for i in [xx]*; do 
        mv $i "$f-${i#*xx}.dump"; done                  && echo "   > Renamed Finished"
  done
  # Return to previous folder
  cd ..
done

#--- Return to previous wd ----
echo 
echo "#============================"
echo 
cd $oldpwd
