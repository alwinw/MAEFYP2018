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

# Copy session files to the ouput location
output=../session-output/
cp -r */ $output
cd $output
echo $(pwd)

#--- Loop Over Airfoils in Folder ----
# Loop over airfoils
for airfoil in */; do
  echo
  echo "#============================"
  # Echo airfoil name
  echo " Found $airfoil"
  # Navigate into folder
  cd $airfoil
  echo " > Entering $(pwd)"
  #--- Loop Over Files in Folder ----
  # Loop over session files
  for sessionfile in *.sesh; do
    echo "#----------------------------"
    echo "   $(date +%D-%T.%10N)"
    # Remove extension
    f="${sessionfile%.*}"                               && echo "   Found $f"               &&
      # Copy to extensionless file
      cp $sessionfile "$f"                              && echo "   > Copy Made"            &&
      # Generate the mesh
      meshpr $f > $f.msh                                && echo "   > Meshpr Finished"      &&
      # Add extra mesh info
      meshpr -i $f > $f.mshi                            && echo "   > Meshpr -i Finished"   &&
      # Enumerate
      enumerate $f > $f.num                             && echo "   > Enumerate Finished"   &&
      # dns (verbosity -v turned off);
      dns $f | grep "Divergence Energy:"                && echo "   > DNS Finished"         &&
      # Add vorticity field t
      addvortfield -G -s $f $f.fld > $f.Gfld            && echo "   > VortGen Added"        &&
      # Convert to ascii
      convert $f.Gfld > $f.flddump                      && echo "   > Convert Finished"     &&
      # Generate wall grad
      wallgrad $f $f.msh > $f.wallgrad                  && echo "   > Wall Grad Finished"   &&
      # Split the output file
      csplit -z $f.flddump /Session/ '{*}' >/dev/null   && echo "   > Split Finished"       &&
      # Rename split files
      for i in [xx]*; do 
        mv $i "$f-${i#*xx}.dump"; done                  && echo "   > Renamed Finished"     &&
      # Remove unneded dump file
      rm $f.msh $f.num $f.fld $f.Gfld $f.mdl $f.flddump $f.his $f.flx
  done
  # Return to previous folder
  cd ..
done

#--- Return to previous wd ----
echo 
echo "#============================"
echo 
cd $oldpwd
