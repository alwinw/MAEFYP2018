# Initial Consultations
 - Various papers to read (and understand)
 - In the stopping phase there is a a pair of LE vortices
 - Vorticity generation (VGen) is due to the pressure gradient an the boundary/wall acceleration
 - VGen could be influenced by the location, time and Re {1,000 - 50,000 approx}
 - Using different reference frames (stationary vs accelerating) should make no difference on the results
 - Melissa's work "validated" the use of Semtex for investigating vorticity
 - Consider investigating a non-rotating cylinder, airfoil at AoA = 0 to compare VGen
 
# Initial Semtex Investigations
 - Try run sp5 (Melissa's file)
 - Try used tecplot macros to create animations
 - VisIt, Paraview are alternatives to Tecplot
 - To generate vorticity use addfield {u, v, w, r, s t}

# Mesh & Interpolation
 - Looked into conformal mapping, however not used since the N-S eqs solved on original domain
 - Airfoil spline created by using the distance between subsequent points

# Gradient Fields
 - Tried mapping to a s-z plane and determining gradients there
 - Tried finite difference methods that requried interpolation

# Semtex Utilities
 - Used wallmesh to extract the wall surface coordinates
 - prob, hist, checkpoint can all be used to extract information

# Semtex Development
 - wallmsh  : add outward normals
 - addfield : add gradients of pressure and vorticity
 - For wallmsh, use sideGeom in order to determine the gradients
 - For addfield, need to add additional fields: klmno

# WallGrad Integration
 - At the corners of two quads, the normals are not exactly equal
 - The average has been taken
 - When the direction of s was changed, so were the directions of the normals in R causing problems
 - Instead of using a join between mesh and field data, try match elements

ctrl+shift+m gives %>% 
