#!/bin/bash

#============================
# Rough Working
# Alwin Wang
#----------------------------

echo "Rough Working"

#--- Working Directories ----
# Save old working dir
oldpwd=$(pwd)
echo $oldpwd
# Navigate to Folder location
folder=../session-files/NACA0012-AoA04/
cd $folder

#--- Loop Over Files in Folder ----
# Session files
sessionfiles=*.sesh
# Loop over session files
for sessionfile in $sessionfiles; do
  # Remove extension
  f="${sessionfile%.*}"                     && echo " Found $f"             &&
    # Copy to extensionless file
    cp $sessionfile "$f"                    && echo "   Copy Made"          &&
    # Generate the mesh
    meshpr $f > "$f.msh"                    && echo "   Meshpr Finished"    &&
    # Enumerate
    enumerate $f > "$f.num"                 && echo "   Enumerate Finished" &&
    # Compare
    compare $f > "$f.rst"                   && echo "   Compare Finished"   &&
    # Generate wall mesh
    wallmesh $f "$f.msh" > "$f.wallmsh"     && echo "   Wall Mesh Finished" &&
    # dns (verbosity -v turned off);
    # use > /dev/null to hide stdout and 2>&1 to hide stderror as well
    dns $f | grep "Divergence Energy:"      && echo "   DNS Finished"       &&
    # Add vorticity field t
    # Note: new field file CANNOT be the same as the original one (else 0)
    addfield -v -s $f "$f.fld" > "$f.vfld"  && echo "   Vorticity Added"    &&
    # Convert to ascii
    convert "$f.vfld" > "$f.dump"           && echo "   Convert Finished"   &&
    # Split the output file
    csplit -z "$f.dump" /Session/ {*}       && echo "   Split Finished"     &&
    # Rename split files
    for i in [xx]*; do 
      mv $i "$f${i#*xx}-.dump"; done         && echo "   Renamed Finished"
done

# Return to previous wd
cd $oldpwd
