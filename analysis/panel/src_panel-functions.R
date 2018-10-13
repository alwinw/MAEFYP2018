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
VLIndVel <- function(panel, panels, output = "Inf") {
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
  # Determine induced velocities in panel coordinates
  indvel <- indvel %>% 
    mutate(
      r1 = sqrt(xcl^2 + ycl^2),
      r2 = sqrt((xcl-x2l)^2 + ycl^2),
      t1 = atan2(ycl, xcl),
      t2 = atan2(ycl, xcl-x2l)) %>% 
    mutate( # Simplified since x1l = 0
      pha = (-ycl/(4*pi*x2l))*((x2l-xcl)*log(r1^2/r2^2)+x2l) - 
        (1/(4*pi*x2l))*((-xcl^2+2*x*x2l+ycl^2)*(t1-t2)+x2l^2*t2),
      phb = (-ycl/(4*pi*x2l))*((    xcl)*log(r1^2/r2^2)-x2l) - 
        (1/(4*pi*x2l))*(( xcl^2        -ycl^2)*(t1-t2)+x2l^2*t2),
      ual = (-ycl*log(r2/r1) + (x2l-xcl)*(t2-t1))      /(2*pi*x2l),
      ubl = ( ycl*log(r2/r1) +       xcl*(t2-t1))      /(2*pi*x2l),
      wal = ((x2l-xcl)*log(r2/r1) - x2l + ycl*(t2-t1)) /(2*pi*x2l),
      wbl = (      xcl*log(r2/r1) + x2l - ycl*(t2-t1)) /(2*pi*x2l) )
  # Transform back into global coordinates
  indvel <- indvel %>% 
    mutate(
      ua = CoRot(ual, wal, 0, 0, theta, "x"),
      wa = CoRot(ual, wal, 0, 0, theta, "y"),
      ub = CoRot(ubl, wbl, 0, 0, theta, "x"),
      wb = CoRot(ubl, wbl, 0, 0, theta, "y") )
  # Influence matrices
  IC <- indvel
  IC[nrow(IC)+1,] = 0 # Added for convinience in lead/lag
  IC <- data.frame(
    u   = IC$ua  + lag(IC$ub,  default = 0),
    w   = IC$wa  + lag(IC$wb,  default = 0),
    phi = IC$pha + lag(IC$phb, default = 0) )
  # Influence row vector
  Kj = IC$u*panel$nx + IC$w*panel$ny
  Lj = IC$u*panel$tx + IC$w*panel$ty
  Pj = IC$phi
  
  # Return output
  if (output == "Inf") 
    return(data.frame(Kj = Kj, Lj = Lj, Pj = Pj))
  else if (output == "ab")
    return(select(indvel, pha, phb, ua, wa, ub, wb))
  else
    return(indvel)
}

# Influence matrices
VLInfMat <- function(panels) {
  Ni = nrow(panels)
  # Apply for each collication point
  panel_li <- split(panels, panels$i)
  # Determine influence matrices
  inf  <- lapply(panel_li, function(panel) VLIndVel(panel, panels))
  Kmat <- t(bind_rows(lapply(inf, function(inf) inf$Kj)))
  Lmat <- t(bind_rows(lapply(inf, function(inf) inf$Lj)))
  Pmat <- t(bind_rows(lapply(inf, function(inf) inf$Pj)))
  # Add kutta condition to Kmat
  Kmat =  rbind(Kmat, c(1, rep(0, Ni-1), 1))
  # Return output
  return(list(Kmat = Kmat, Lmat = Lmat, Pmat = Pmat))
}

#--- Panel Method Solution                                        ----
# aoa = 4
# U = 1
VLSol <- function(coord, aoa, U) {
  # Determine influence matrices
  panels <- VLPanel(coord)
  infmat <- VLInfMat(panels)
  Kmat   <- infmat$Kmat
  Lmat   <- infmat$Lmat
  Pmat   <- infmat$Pmat
  # Determine RHS
  alpha = aoa*pi/180 
  Uinf = U*cos(alpha)
  Vinf = U*sin(alpha)
  bvec = -(Uinf*panels$nx + Vinf*panels$ny)
  bvec = c(bvec, 0)
  # Solve for gamma
  gvec = as.matrix(solve(Kmat, bvec))
  # Induced velocities
  vvec = Lmat %*% gvec + Uinf*panels$tx + Vinf*panels$ty
  cp   = 1 - vvec^2
  
  return(list(g = gvec, v = vvec, cp = cp))
}

