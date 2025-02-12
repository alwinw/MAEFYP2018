# Create required setup files
# Run dns on session file
mkdir results
mkdir images
for sessionfile in *.sesh; do
  f="${sessionfile%.*}"                   &&
  cp $sessionfile  results/$sessionfile   &&
  cp bndry_prf.dat results/bndry_prf.dat  &&
  cd results                              &&
  cp $sessionfile  $f                     &&
  meshpr    $f > $f.msh                   &&
  meshpr -i $f > $f.mshi                  &&
  enumerate $f > $f.num                   &&
  dns $f | grep "Divergence Energy:"      &&
  addvortfield -G -s $f $f.fld > $f.Gfld  &&
  convert  $f.Gfld   > $f.flddump         &&
  compare  $f $f.fld > /dev/null          &&
  wallmesh $f $f.msh > $f.wallmsh         &&
  wallgrad $f $f.msh > $f.wallgrad        &&
  csplit -z "$f.flddump" /Session/ '{*}' >/dev/null &&
  for i in [xx]*; do 
    mv $i "$f-${i#*xx}.dump"; done                  &&
  echo $sessionfile complete
  cd ..
done

