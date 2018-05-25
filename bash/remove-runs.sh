# /bin/bash

#============================
# Remove Runs from Rep
# Alwin Wang
#----------------------------
# https://stackoverflow.com/questions/48699293/how-to-i-disable-git-lfs
# https://github.com/git-lfs/git-lfs/issues/910#issuecomment-238389388

# git lfs ls-files | cut -d ' ' -f 3 > lfs-files.txt
# cat bash/lfs-files.txt | xargs touch

# git rm -r --cached .
# git add .

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
      # Remove unnecessary files
      rm "$f"                                           && echo "   > Removed session"    &&
      # Remove the mesh
      rm "$f.msh"                                       && echo "   > Removed Meshpr"     &&
      # Remove extra mesh info
      rm "$f.mshi"                                      && echo "   > Removed Meshpr -i"  &&
      # Remove Enumerate
      rm "$f.num"                                       && echo "   > Removed Enumerate"  &&
      # Remove Compare
      rm "$f.rst"                                       && echo "   > Removed Compare"    &&
      # Remove wall mesh
      rm "$f.wallmsh"                                   && echo "   > Removed Wall Mesh"  &&
      # Remove vorticity field
      rm "$f.vfld"                                      && echo "   > Vorticity Added"    &&
      # Remove Convert to ascii
      rm "$f.flddump"                                   && echo "   > Remove Convert"     &&
      # Remove split files
      rm "*.dump"                                       && echo "   > Remove Split"
  done
  # Return to previous folder
  cd ..
done

#--- Return to previous wd ----
echo 
echo "#============================"
echo 
cd $oldpwd
