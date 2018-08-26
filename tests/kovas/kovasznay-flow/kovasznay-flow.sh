# Kovasznay Flow

meshpr kovas1 > kovas1.msh
enumerate kovas1 > kovas1.num
compare kovas1 > kovas1.rst
dns kovas1
addfield     -v -s kovas1 kovas1.fld > kovas1.vfld
addvortfield -G -s kovas1 kovas1.fld > kovas1.Gfld
convert kovas1.vfld > kovas1.vflddump
convert kovas1.Gfld > kovas1.Gflddump
compare kovas1 kovas1.fld > /dev/null
