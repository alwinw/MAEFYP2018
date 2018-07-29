# wallmesh.cpp      | 2018.05.06
 - Changed wall mesh to double precision (patch by Hugh B)

# meshpr.cpp        | 2018.05.06
 - Added `-i` to options to print extra information

# wallgrad.cpp      | 2018.05.29
 - Created wallgrad.cpp to calculate gradients at wall

    grep -r "sideGeom" "."

 - /src/element.cpp and /src/edge.cpp

# addvortfield.cpp  | 2018.07.23
 - Added additional auxfields for gradients
 - `-G` for VORTGEN option
 - Debugging using gdb
 
    $ gdb addvortfield
    (gdb) run -G -s kovas1 kovas1.fld
    (gdb) backtrace
 
 - Example
 
    addvortfield -G -s kovas1 kovas1.fld > kovas.vfld
    