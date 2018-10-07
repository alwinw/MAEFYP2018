#============================#
# Panel Method
# Alwin Wang
#----------------------------#
library(ggplot2)
library(dplyr)
library(tidyr)
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 

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

VLInf <- function(panel, panels) {
  # Need to rotate the collocation point into the panel coordinate
  indvel <- panels %>% 
    mutate(
      xc = panel$x,
      yc = panel$y) %>% 
    mutate(
      x1d = CoRot(x1, y1, x1, y1, -theta, "x"),
      y1d = CoRot(x1, y1, x1, y1, -theta, "y"),
      x2d = CoRot(x2, y2, x1, y1, -theta, "x"),
      y2d = CoRot(x2, y2, x1, y1, -theta, "y"),
      xcd = CoRot(xc, yc, x1, y1, -theta, "x"),
      ycd = CoRot(xc, yc, x1, y1, -theta, "y"))
  # ggplot(indvel) +
  #   geom_segment(aes(x1d, y1d, xend = x2d, yend = y2d)) +
  #   geom_point(aes(xcd, ycd, colour = i))
  # 
  # Determined induced velocities in panel coordinates
  indvel <- indvel %>% 
    mutate(
      r1 = sqrt((xcd-x1d)^2 + ycd^2),
      r2 = sqrt((xcd-x2d)^2 + ycd^2),
      t1 = atan2(ycd, xcd-x1d),
      t2 = atan2(ycd, xcd-x2d)) %>% 
    mutate(
      t1 = ifelse(t1<0, t1+2*pi, t1),
      t2 = ifelse(t2<0, t2+2*pi, t2) ) %>% 
    mutate(
      ua = ycd/(2*pi)*(-1/(x2d-x1d))*log(r2/r1) + ((x2d-x1d)-1*(xcd-x1d))/(2*pi*(x2d-x1d))*(t2-t1),
      wa = ((x2d-x1d)-1*(xcd-x1d))/(2*pi*(x2d-x1d))*log(r2/r1) + 1/(2*pi)*(-1/(x2d-x1d))*((x2d-x1d)+(t1-t2)*ycd) ) %>% 
    mutate(
      ub = ycd/(2*pi)*(+1/(x2d-x1d))*log(r2/r1) + (1*(xcd-x1d))/(2*pi*(x2d-x1d))*(t2-t1),
      wb = (1*(xcd-x1d))/(2*pi*(x2d-x1d))*log(r2/r1) + 1/(2*pi)*(1/(x2d-x1d))*((x2d-x1d)+(t1-t2)*ycd) )
  
  # Transform back into global coordinates
  indvel <- indvel %>% 
    mutate(
      ua = CoRot(ua, wa, 0, 0, theta, "x"),
      wa = CoRot(ua, wa, 0, 0, theta, "y"),
      ub = CoRot(ub, wb, 0, 0, theta, "x"),
      wb = CoRot(ub, wb, 0, 0, theta, "y") ) %>% 
    select(ua, wa, ub, wb)
  indvel[nrow(indvel)+1,] = 0 # Added for convinience in lead/lag
  
  indvel <- data.frame(
    u = indvel$ua + lag(indvel$ub, default = 0),
    w = indvel$wa + lag(indvel$wb, default = 0) )
  
  # Influence row
  Kj = indvel$u*panel$nx + indvel$w*panel$ny
  
  return(Kj)
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
  
  # gamma = c(-0.8862, -0.9905, 0.4848, 1.233, 0.8862)
  
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

coord <- data.frame(
  x = c(1,  0.3,  0,  0.3, 1),
  y = c(0, -0.05, 0, 0.05, 0))
panels <- VLPan(coord)
vlsoln <- VLSol(panels)
ggplot(coord, aes(x, y)) + geom_path() + geom_point() + 
  coord_fixed() +
  geom_segment(aes(x, y, xend = x+nx, yend = y+ny), panels)



# Chord Vector
cvec = seq(-1, 1, length.out = 10)
cvec = AirfoilSamp(cvec, 1.5)
coord = AirfoilCoord(cvec, 4412)
ggplot(coord$coord, aes(x, y)) + geom_path() + geom_point() + 
  geom_path(aes(x, y), coord$camb) +
  coord_fixed()
max(coord$coord$y) - min(coord$coord$y)
panels <- VLPan(coord$coord)
vlsoln <- VLSol(panels)
ggplot(vlsoln, aes(x, -Cp)) + geom_path() + ylim(-2, 1)

