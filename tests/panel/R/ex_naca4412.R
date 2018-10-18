#============================#
# NACA4412 Example
# Alwin Wang
#----------------------------#

library(dplyr)
library(tidyr)

#--- Transform a set of (x,y) based on an AoA (in deg) ----
AoATransform <- function(data, AoA) {
  # Assumes X is column 1 and Y is column 2
  # Store the input data
  odata <- data
  ocolnames <- colnames(data)
  colnames(data[c(1,2)]) <- c("x", "y")
  # Apply the transformation
  data <- select(data, x, y) %>%
    mutate(
      r = sqrt(x^2 + y^2),
      theta = atan(y/x),
      theta = ifelse(is.na(theta),0,theta),
      theta = ifelse(x < 0, theta + pi, theta),
      theta = theta - AoA*pi/180,
      x = r*cos(theta),
      y = r*sin(theta)
    ) %>%
    select(x, y)
  # Re-combine the new data with the old data and restore colnames
  odata[c(1,2)] = data
  colnames(odata) <- ocolnames
  return(odata)
}

#--- Surface Coordinates for a NACA 4 digit airfoil ----
AirfoilCurve <- function(x = 0, out = "all") {
  # Test if x is within range
  on = ifelse(x >= a & x <= (a + c), 1, 0) # allows for root-finding
  # Determine the camber line yc
  yc = ifelse(x < p * c + a, 
              m/p^2 * (2*p*((x-a)/c) - ((x-a)/c)^2),
              m/(1-p)^2  * (1 - 2*p + 2*p*((x-a)/c) - ((x-a)/c)^2)
  )
  # Determine the gradient of the camber line dycdx
  dycdx = ifelse(x < p * c + a,
                 2*m/p^2 * (p - (x-a)/c),
                 2*m/(1-p)^2 * (p - (x-a)/c)
  )
  # Determine the magnitude and direction of the thickness
  theta = atan(dycdx)
  yt = 5*t*(0.2969*sqrt(abs((x-a)/c)) - 0.1260*((x-a)/c) - 0.3516*((x-a)/c)^2 +
              0.2843*((x-a)/c)^3 - 0.1036*((x-a)/c)^4)
  # Add the thickness to the camber line
  xU = x - yt*sin(theta)
  yU = yc + yt*cos(theta)
  xL = x + yt*sin(theta)
  yL = yc - yt*cos(theta)
  # Output depending on the Out parameter
  if(out == "all")
    summary = data.frame(x, yc, dycdx, theta, yt,  xU, yU,  xL, yL)
  else if(out == "coord")
    summary = data.frame(x, xU, yU, xL, yL)
  else if (out == "upper")
    summary = data.frame(x = xU, y = yU)
  else if (out == "lower")
    summary = data.frame(x = xL, y = yL)
  # Return the output
  return(summary * on)
}

#--- Outputs Airfoil Points into (x,y) columns and AoA transform for plotting ----
AirfoilCoord <- function(xmin = a, xmax = c + a, AoA = 0, res = 100) {
  # Cluster points around LE and TE
  xvec = abs(a) * sin(seq(xmin, xmax, length.out = res)*pi/c)
  # Generate coordinates in a tidy format
  coord = AirfoilCurve(xvec, out = "coord") %>%
    rename(xO = x) %>%
    gather(key, value, -xO) %>%
    mutate(coord = substr(key,1,1), surf = substr(key, 2,2)) %>%
    select(-key) %>%
    spread(coord, value) %>%
    mutate(surf = factor(surf, levels = c("U", "L"))) %>%
    arrange(surf, xO*ifelse(surf=="U", 1, -1)) %>%
    select(x, y, surf)
  coord = AoATransform(coord, AoA = AoA)
  return(coord)
}

#--- Function for better sampling of points ----
AirfoilSamp <- function(xvec, del = c*8e-6, cylinder = FALSE) {
  # xvec = seq(a, a+c, by = 0.01)
  # Sample according to a cubic function
  xvec = -2*a/c^3 * (xvec - a)^3 + a
  # Add extra x values for interpolation
  if (cylinder != FALSE & xvec[1] == a) {
    # Determine the number of points from -theta_c to theta_c
    xadd = seq(-0.0001, -thetac,  
               length.out = ceiling(length(xvec[xvec < xsamp])/2 + 1))
    xadd = xadd[xadd != -thetac]
    # 'encode it' and combine
    xadd = a - abs(a) + xadd
    # Return the result depending on what's required
    if (cylinder == TRUE)
      xvec = c(xadd, xvec)
    if (cylinder == "only")
      xvec = xadd
  }
  
  xvec = xvec[xvec != a]
  
  # Adjust the TE value
  if (xvec[length(xvec)] == a + c)
    xvec[length(xvec)] = a + c - sign(a + c)*abs(a + c)*del
  
  return(xvec)
}



NACA = 4412
# Max camber; Location of m; Thickness
m = (NACA %/% 1000) / 100
p = (NACA %/% 100 %% 10) / 10
t = (NACA %% 100) / 100
# Chord; x-shift
c = 1
a = - 1/2

AirfoilSamp(seq(a, a+c, by = 0.01))

AirfoilCoord(a, c+a)
