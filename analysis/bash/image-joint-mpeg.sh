#!/bin/bash

#============================
# Loop over folder outputs to create animations
# Alwin Wang
#----------------------------

# Save old working dir
oldpwd=$(pwd)
echo $oldpwd
# Navigate to Folder location
folder=../plot-output/
cd $folder
echo $(pwd)

#--- Loop Over Cases in Folder ----
# Create a named pipe
for case in */; do
  echo
  echo "  Found $case"
  cd $case
  for image in *_ns.png; do
    # Remove trailing information
    f="${image%_*}"
    echo "  "$f
    /usr/bin/convert \( $f"_ns.png" $f"_err.png" \) -append $f"_nscom.png" &
    /usr/bin/convert \( $f"_le.png" $f"_a.png" \)   -append $f"_lecom.png" &
    wait
  done


  # N.B. case includes a trailing slash, remove it using ${..%/}
  #/usr/bin/convert -delay 10 -loop 0 *_a.png   m2v:../${case%/}"_a.mpeg" &
  #/usr/bin/convert -delay 10 -loop 0 *_err.png m2v:../${case%/}"_err.mpeg" &
  #wait
  #/usr/bin/convert -delay 10 -loop 0 *_le.png  m2v:../${case%/}"_le.mpeg" &
  #/usr/bin/convert -delay 10 -loop 0 *_ns.png  m2v:../${case%/}"_ns.mpeg" &
  #wait
  cd ..
done

cd $oldpwd
