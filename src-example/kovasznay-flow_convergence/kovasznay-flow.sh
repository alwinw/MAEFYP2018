# Create required setup files
# Create session files 
N_P=( 3 4 5 6 7 8 9 10 11 12 13 14 15 )
# echo $N_P[*]
for i in $N_P[*]; do
  name=$(printf "%02d" i)
  echo N_P = $i
  sed  "s/N_P =/N_P    = $i/g" kovas >> results/kovas$name.sesh
done
cd results
# Run dns on them
for sessionfile in *.sesh; do
  f="${sessionfile%.*}"   &&
  mv $sessionfile $f      &&
  meshpr -i $f > $f.mshi  &&
  enumerate $f > $f.num   &&
  dns $f | grep "Divergence Energy:"      &&
  addvortfield -G -s $f $f.fld > $f.Gfld  &&
  convert $f.Gfld > $f.flddump            &&
  compare $f $f.fld > /dev/null           &&
  rm $f.fld $f.Gfld $f.mdl $f.flx $f.num  &&
  echo $sessionfile complete
done
cd ..
