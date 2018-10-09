#============================#
# Panel Method
# Alwin Wang
#----------------------------#
library(ggplot2)
library(dplyr)
library(tidyr)
theme_set(theme_bw())                                           # Set black and white theme

AirfoilSamp <- function(cvec, p) {
  a = min(cvec)
  c = cvec[length(cvec)] - cvec[1]
  xvec = ((0.5*sin((cvec - 0.5*c - a)*pi/c) + 0.5)^p)*c + a # Asymmetric
  return(xvec)
}

AirfoilCoord <- function(
  cvec, NACA = NULL,
  mc = NULL, p = NULL, tc = NULL) {
  # Determine airfoil properties
  a = min(cvec)
  c = cvec[length(cvec)] - cvec[1]
  if (is.numeric(NACA)) {
    mc = (NACA %/% 1000) / 100
    p = (NACA %/% 100 %% 10) / 10
    tc = (NACA %% 100) / 100
  }
  t = tc*c
  m = mc*c
  # Generate coordinates
  x  = cvec
  yc = ifelse(x < p*c+a, 
              m/p^2*(2*p*((x-a)/c)-((x-a)/c)^2),
              m/(1-p)^2*(1-2*p + 2*p*((x-a)/c)-((x-a)/c)^2))
  dycdx = ifelse(x < p*c+a,
                 2*m/p^2*(p-(x-a)/c),
                 2*m/(1-p)^2*(p-(x-a)/c))
  theta = atan(dycdx)
  yt = 5*t*(0.2969*sqrt(abs((x-a)/c)) - 0.1260*((x-a)/c) - 0.3516*((x-a)/c)^2 +
              0.2843*((x-a)/c)^3 - 0.1036*((x-a)/c)^4)
  xU = x - yt*sin(theta)
  yU = yc + yt*cos(theta)
  xL = x + yt*sin(theta)
  yL = yc - yt*cos(theta)
  # Output
  coord <- data.frame(x = c(rev(xL), xU), y = c(rev(yL), yU)) %>% 
    unique(.)
  rownames(coord) <- NULL
  return(list(coord = coord, camb = data.frame(x = x, y = yc)))
}

VLPan <- function(coord) {
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
           theta = atan2(dy, dx)) %>% 
    select(-nl)
  # Return the output
  return(panels)
}

CoRot <- function(x, y, x0, y0, theta, var) {
  if (var == "x") 
    return(cos(theta)*(x-x0) - sin(theta)*(y-y0))
  if (var == "y") 
    return(sin(theta)*(x-x0) + cos(theta)*(y-y0))
}

VLInf <- function(panel, panels, vel = FALSE) {
  # Need to rotate the collocation point into the panel coordinate
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
      ycl = CoRot(xc, yc, x1, y1, -theta, "y"))

  # Determined induced velocities in panel coordinates
  indvel <- indvel %>% 
    mutate(
      r1 = sqrt(xcl^2 + ycl^2),
      r2 = sqrt((xcl-x2l)^2 + ycl^2),
      t1 = atan2(ycl, xcl),
      t2 = atan2(ycl, xcl-x2l)) %>% 
    mutate(
      t1 = ifelse(t1<0, t1+2*pi, t1),
      t2 = ifelse(t2<0, t2+2*pi, t2) ) %>%
    mutate( # Simplified since x1l = 0
      ual = (-ycl*log(r2/r1) + (x2l-xcl)*(t2-t1))      /(2*pi*x2l),
      ubl = ( ycl*log(r2/r1) +       xcl*(t2-t1))      /(2*pi*x2l),
      wal = ((x2l-xcl)*log(r2/r1) - x2l + ycl*(t2-t1)) /(2*pi*x2l),
      wbl = (      xcl*log(r2/r1) + x2l - ycl*(t2-t1)) /(2*pi*x2l) )
  
  # Transform back into global coordinates
  indvel <- indvel %>% 
    # mutate(
    #   ua = CoRot(ual, wal, 0, 0, theta, "x"),
    #   wa = CoRot(ual, wal, 0, 0, theta, "y"),
    #   ub = CoRot(ubl, wbl, 0, 0, theta, "x"),
    #   wb = CoRot(ubl, wbl, 0, 0, theta, "y") ) %>% 
    # select(ua, wa, ub, wb)
    select(ual, wal, ubl, wbl, t1, t2, ycl)
  
  # IC <- indvel
  # IC[nrow(IC)+1,] = 0 # Added for convinience in lead/lag
  # 
  # IC <- data.frame(
  #   u = IC$ua + lag(IC$ub, default = 0),
  #   w = IC$wa + lag(IC$wb, default = 0) )
  # 
  # # Influence row
  # Kj = IC$u*panel$nx + IC$w*panel$ny
  
  if (vel) return(indvel)
  else return(Kj)
}

VLSol <- function(panels) {
  Ni = nrow(panels)
  # Apply for each collication point
  panel_li <- split(panels, panels$i)
  K <- lapply(panel_li, function(panel) VLInf(panel, panels))
  K =  t(bind_rows(K))
  K =  rbind(K, c(1, rep(0, Ni-1), 1))
  # RHS
  alpha = 5*pi/180
  Uinf = cos(alpha)
  Vinf = sin(alpha)
  b =  -(Uinf*panels$nx + Vinf*panels$ny)
  b =  c(b, 0)
  # Solve
  gamma = solve(K, b)
  
  # Output properties
  # CHORD = left most and right most EucDist
  panels$ga = gamma[1:Ni]
  panels$gb = gamma[2:(Ni+1)]
  panels <- panels %>%
    mutate(U = Uinf*dx + Vinf*dy + (ga+gb)/4) %>% 
    mutate(Cp = 1-(U/Uinf)^2) %>% 
    mutate(Cl = (ga+gb)*c/1)
  
  return(panels)
}

# Diamond
coord <- data.frame(
  x = c(1,  0.3,  0,  0.3, 1),
  y = c(0, -0.05, 0, 0.05, 0))
panels <- VLPan(coord)
vlsoln <- VLSol(panels)
ggplot(coord, aes(x, y)) + geom_path() + geom_point() + 
  coord_fixed() +
  geom_segment(aes(x, y, xend = x+nx, yend = y+ny), panels)

# NACA0012.10
# Note: manually added the last point to close the loop
coord <- data.frame(
  x = c(1.000000, 0.904508, 0.654508, 0.345491, 0.095491,
        0.000000, 0.095491, 0.345491, 0.654508, 0.904508,
        1.000000),
  y = c(0.000000,-0.013914,-0.040917,-0.059575,-0.046049,
        0.000000, 0.046049, 0.059575, 0.040917, 0.013914,
        0.000000)
)
panels <- VLPan(coord)
ggplot(coord, aes(x, y)) + geom_path() + geom_point() + 
  coord_fixed() +
  geom_segment(aes(x, y, xend = x+nx, yend = y+ny), panels)
panel_li <- split(panels, panels$i)
options(digits = 5, scipen = -4)
lapply(panel_li, function(panel) VLInf(panel, panels, TRUE))


# Chord Vector
cvec = seq(-1, 1, length.out = 100)
cvec = AirfoilSamp(cvec, 1.5)
coord = AirfoilCoord(cvec, 4412)
ggplot(coord$coord, aes(x, y)) + geom_path() + geom_point() + 
  geom_path(aes(x, y), coord$camb) +
  coord_fixed()
max(coord$coord$y) - min(coord$coord$y)
panels <- VLPan(coord$coord)
vlsoln <- VLSol(panels)
ggplot(vlsoln, aes(x, Cp)) + geom_path() + ylim(-2, 1) + scale_y_reverse()

