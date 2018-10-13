#============================#
# Panel Method Functions
# Alwin Wang
#----------------------------#

# Panel method: Vortex Panels of Linearly Varying Strength
# Kuethe & Chow: 
# * Section 5.10 The Airfoil of Arbitrary Thickness and Camber
# > eq 5.43 for velocity potential
# > pg 161 example code and output
# Katz & Plotkin: 
# * Section 10.3.3 LInear Vortex Distribution
# > eq 10.71 for velocity potential
# > eq 10.72, 10.73 for u and w
# * Section 11.4.2 Linear-Strength Vortex Method

#--- Notation                                                     ----
# i       = (x,y)_i, collocation point
# j       = gamma_j, linear strength vortex
# Ni      = number of collocation points (and panels)
# Nj      = number of gamma
# (x,y)   = collocation point
# (x1,y1) = 'left' panel end
# (x2,y2) = 'right' panel end
# L       = local panel coordinate system
# G       = global coordinate system

#--- Panel Geometry                                               ----
# Determine the geometry of the panels
VLPanel <- function(coord) {
  # Num of collocation points
  Ni = nrow(coord) - 1
  # Determine panel segments
  panels <-  data.frame(
    i  = 1:(Ni+1),
    x  = (lead(coord$x) + coord$x)/2,
    y  = (lead(coord$y) + coord$y)/2,
    x1 = coord$x,
    y1 = coord$y,
    x2 = lead(coord$x),
    y2 = lead(coord$y)) %>% 
    mutate (
      c  = sqrt((x2-x1)^2 + (y2-y1)^2),
      dx = x2 - x1,
      dy = y2 - y1)
  panels   <- panels[1:Ni,]
  # Determine normal components using cross product 
  # Clockwise: (0, 0, -1) | AntiClock: (0, 0, 1) 
  panels <- panels %>% 
    mutate(nx = -dy, ny = dx,
           nl = sqrt(nx^2+ny^2),
           nx = nx/nl,
           ny = ny/nl,
           tx = dx/c,
           ty = dy/c,
           theta = atan2(dy, dx)) %>% 
    select(-nl)
  # Return the output
  return(panels)
}

#--- Influence Matrices                                           ----
# Coordinate transform between local and global coordinate systems
CoRot <- function(x, y, x0, y0, theta, var) {
  if (var == "x") 
    return(cos(theta)*(x-x0) - sin(theta)*(y-y0))
  if (var == "y") 
    return(sin(theta)*(x-x0) + cos(theta)*(y-y0))
}

# Induced velocity and potentials
VLIndVel <- function(panel, panels) {
  # Transform the panels into local coordinates of panel
  indvel <- panels %>% 
    mutate(
      xc = panel$x,
      yc = panel$y) %>% 
    mutate( # Due to coordinate transform x1l = y1l = y2l = 0
      x1l = 0,
      y1l = 0,
      x2l = CoRot(x2, y2, x1, y1, -theta, "x"),
      y2l = 0, 
      xcl = CoRot(xc, yc, x1, y1, -theta, "x"),
      ycl = CoRot(xc, yc, x1, y1, -theta, "y")) %>% 
    mutate( # Remove floating point error of atan(0)
      ycl = ifelse(i == panel$i, 0, ycl))
  # Determine induced velocities
  # Determined induced velocities in panel coordinates
  indvel <- indvel %>% 
    mutate(
      r1 = sqrt(xcl^2 + ycl^2),
      r2 = sqrt((xcl-x2l)^2 + ycl^2),
      t1 = atan2(ycl, xcl),
      t2 = atan2(ycl, xcl-x2l)) %>% 
  
}