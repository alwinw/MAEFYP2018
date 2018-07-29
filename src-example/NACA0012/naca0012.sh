# Create required setup files
# Run dns on session file
for sessionfile in *.sesh; do
  f="${sessionfile%.*}"                   &&
  cp $sessionfile  results/$f             &&
  cp bndry_prf.dat results/bndry_prf.dat  &&
  cd results                              &&
  meshpr    $f > $f.msh                   &&
  meshpr -i $f > $f.mshi                  &&
  enumerate $f > $f.num                   &&
  dns $f | grep "Divergence Energy:"      &&
  addvortfield -G -s $f $f.fld > $f.Gfld  &&
  convert  $f.Gfld   > $f.flddump         &&
  compare  $f $f.fld > /dev/null          &&
  wallgrad $f $f.msh > $f.wallgrad        &&
  csplit -z "$f.flddump" /Session/ '{*}' >/dev/null &&
  for i in [xx]*; do 
    mv $i "$f-${i#*xx}.dump"; done                  &&
  echo $sessionfile complete
  cd ..
done

