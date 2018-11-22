# Create required setup files
# Run dns on session file
mkdir -p results
mkdir -p images
mkdir -p tecplot
for sessionfile in *.sesh; do
  f="${sessionfile%.*}"                   &&
  cp $sessionfile  results/$sessionfile   &&
  cd results                              &&
  cp $sessionfile  $f                     &&
  meshpr    $f > $f.msh                   &&
  meshpr -i $f > $f.mshi                  &&
  massmat   $f > $f.mass                  &&
  enumerate $f > $f.num                   &&
  dns $f | grep "Divergence Energy:"      &&
  addvortfield -G -s $f $f.fld > $f.Gfld  &&
  convert  $f.Gfld   > $f.flddump         &&
  compare  $f $f.fld > /dev/null          &&
  wallmesh $f $f.msh > $f.wallmsh         &&
  wallgrad $f $f.msh > $f.wallgrad        &&
  csplit -z "$f.flddump" /Session/ '{*}' >/dev/null            &&
  for i in [xx]*; do 
    dump="$f-${i#*xx}"                    &&
    mv $i $dump.dump                      &&
    time=$(sed -n '5p' $dump.dump | grep -Eo "[0-9]+\.[0-9]+") &&
    echo $time                            &&
    sem2tec -n 0 -o $dump.dat -m $f.msh $dump.dump             &&
    sed -i "/POINT/s/$/, SOLUTIONTIME=$time/" $dump.dat        &&
    preplot $dump.dat > /dev/null         &&
    rm $dump.dat                          &&
    mv $dump.plt ../tecplot/$dump.plt     &&
  done                                    &&
  echo $sessionfile complete
  cd ..
done

