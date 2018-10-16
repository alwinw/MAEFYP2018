# Create required folders (if necessary)
mkdir -p results
mkdir -p images
# Create required setup files
# Run dns on session file
for sessionfile in *.sesh; do
  f="${sessionfile%.*}"                   &&
  cp $sessionfile  results/$sessionfile   &&
  cp bndry_prf.dat results/bndry_prf.dat  &&
  cd results                              &&
  cp $sessionfile  $f                     &&
  meshpr    $f > $f.msh                   &&
  meshpr -i $f > $f.mshi                  &&
  massmat   $f > $f.mass                  &&
  enumerate $f > $f.num                   &&
  dns $f | grep "Divergence Energy:"      &&
  addvortfield -G -s $f $f.fld > $f.Gfld  &&
  convert  $f.Gfld    > $f.flddump        &&
  compare  $f $f.fld  > /dev/null         &&
  wallmesh $f $f.msh  > $f.wallmsh        &&
  wallgrad $f $f.msh  > $f.wallgrad       &&
  csplit -z "$f.flddump" /Session/ '{*}' >/dev/null   &&
  for i in [xx]*; do 
    mv $i "$f-${i#*xx}.dump"; done                    &&
  integral $f $f.Gfld > $f.integral       &&
  csplit -z "$f.integral" /timestep/ '{*}' >/dev/null &&
    for i in [xx]*; do 
    mv $i "$f-${i#*xx}.integral"; done                &&
  echo $sessionfile complete
  cd ..
done

