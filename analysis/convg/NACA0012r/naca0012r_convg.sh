# Create required folders (if necessary)
mkdir -p results
mkdir -p images
# Create array of N_P to loop over
N_P=( 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 )
# N_P=( 3 4 )
for i in $N_P[*]; do
  # Set up dns case
  name=$(printf "%02d" i)
  echo N_P = $name
  mkdir -p results/NACA0012r-N_P$name
  sed "s/N_P = X/N_P       = $i/g" \
    naca0012r >    results/NACA0012r-N_P$name/naca0012r-N_P$name.sesh
  cp bndry_prf.dat results/NACA0012r-N_P$name/bndry_prf.dat
  # Run dns on session file
  cd results/NACA0012r-N_P$name
  for sessionfile in *.sesh; do
    f="${sessionfile%.*}"                   &&
    cp $sessionfile  $f                     &&
    meshpr    $f > $f.msh                   &&
    meshpr -i $f > $f.mshi                  &&
    massmat   $f > $f.mass                  &&
    enumerate $f > $f.num                   &&
    dns $f | grep "Divergence Energy:"      &&
    addvortfield -G -s $f $f.fld > $f.Gfld  &&
    convert  $f.Gfld   > $f.flddump         &&
  # compare  $f $f.fld > /dev/null          &&
  # wallmesh $f $f.msh > $f.wallmsh         &&
    wallgrad $f $f.msh > $f.wallgrad        &&
    csplit -z "$f.flddump" /Session/ '{*}' >/dev/null &&
    for i in [xx]*; do 
      mv $i "$f-${i#*xx}.dump"; done                  &&
    echo $sessionfile complete
    rm $f.msh $f.num $f.fld $f.Gfld $f.mdl $f.flddump $f.his
  done
  cd ../..
done

