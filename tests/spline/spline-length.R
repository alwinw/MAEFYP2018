# Libraries
library("pracma")
library("dplyr")
library("ggplot2")
library("pbapply")

# Functions
LagDist <- function(x, y) sqrt((x - lag(x))^2 + (y - lag(y))^2)
EucDist <- function(x, y) sqrt(x^2 + y^2)

# Spline Length
CubicSplineArcLength <- function(tvec, dcsx, dcsy, method = "Kronrod") {
  # s = int sqrt(dx/dt^2 + dy/dt^2) dt
  splinefun <- function(t) sqrt(ppval(dcsx, t)^2 + ppval(dcsy, t)^2)
  if (method %in% c("Kronrod", "Clenshaw","Simpson")) {
    integral(splinefun, tvec[1], tvec[2], method)
  } else if (method == "Quadgr") {
    quadgr(splinefun, tvec[1], tvec[2])$value
  } else if (method == "Romberg") {
    romberg(splinefun, tvec[1], tvec[2])$value
  } else {
    integrate(splinefun, tvec[1], tvec[2])$value
  }
}

# Determine derivatives and antiderivatives of cubic splines
CubicSplineCalc <- function(cs, order = 0) {
  # Object to return
  ccs <- cs
  # nth Derivative or nth Integral
  while(order != 0) {
    if(order < 0) {
      # Derivative
      if (ncol(ccs$coefs) > 2) ccs$coefs = t(apply(ccs$coefs, 1, polyder))
      else if (ncol(ccs$coefs) == 2) ccs$coefs = apply(ccs$coefs, 1, polyder)
      else ccs$coefs = rep(0, length(ccs$coefs))
      # Ensure correct data structure
      ccs$coefs <- as.matrix(ccs$coefs)
      # Increase order
      order = order + 1}
    else {
      # Integrate
      ccs$coefs = t(apply(ccs$coefs, 1, polyint))
      # Set correct intercepts
      int = c(0)
      for (i in 2:(nrow(ccs$coefs))) {
        int[i] = polyval(ccs$coefs[i-1,], ccs$breaks[i] - ccs$breaks[i-1])
      }
      int = cumsum(int)
      ccs$coefs[,ncol(ccs$coefs)] = int
      # Decrease order
      order = order -1}
  }
  # Adjust to correct order
  ccs$order = ncol(ccs$coefs)
  return(ccs)
}

SplineDist <- function(npoints, method, seed) {
  # Generate random data
  # Set the seed value
  ptm <- proc.time()
  theta = seq(0, 2*pi, length.out = npoints) + runif(50, -0.1, 0.1)
  theta = sort(theta)
  theta =  c(0, theta[theta >=0 & theta <= 2*pi], 2*pi)
  data = data.frame(x = cos(theta), y = sin(theta))
  ggplot(data, aes(x, y)) + geom_path() + geom_point() + coord_fixed()
  # Determine Spline Distance
  data$s = c(0, LagDist(data$x, data$y)[-1])
  data$s = cumsum(data$s)
  print(max(data$s)/(2*pi))
  error = 1                                                       # Initialise variable
  # print("Starting Loop")
  while (abs(error) > 0.001) {                                    # Loop for required error
    csx <- cubicspline(data$s, data$x)                            # Cublic spline
    csy <- cubicspline(data$s, data$y)
    dcsx = CubicSplineCalc(csx, -1)                               # Derivative of cubic splines
    dcsy = CubicSplineCalc(csy, -1)
    cat("Determining spline distance\n")
    s <- cbind(data$s, lead(data$s))[-nrow(data),]                # Interval (s1, s2)
    s <- pbapply(s, 1, CubicSplineArcLength, dcsx, dcsy, method)  # Arc length of interval
    s <- c(0, cumsum(s))                                          # Cumulative arc lenth
    error  = sum(data$s - s)                                      # Error
    data$s = s                                                    # Update loop
    # print(error)
    # print(max(s)/(2*pi))
  }
  ptm <- proc.time() - ptm
  output <- data.frame(
    npoints = npoints,
    method = method,
    s = max(s),
    itererror = error,
    perierror = max(s)/(2*pi),
    timeuser = ptm[1],
    timesysm = ptm[2],
    timeelas = ptm[3]
  )
  return(output)
}

# SplineDist(50, "Kronrod")
# SplineDist(50, "Clenshaw")
# SplineDist(50, "Simpson")
# SplineDist(50, "Quadgr")
# SplineDist(50, "Romberg")
# SplineDist(50, "Base")

npoints <- seq(50, 300, 50)
methods <- c("Kronrod", "Clenshaw", "Simpson", "Quadgr", "Romberg", "Base")

compare <- data.frame(
  npoint = rep(npoints, each = length(methods)),
  method = rep(methods, times = length(npoints)),
  stringsAsFactors = FALSE
)

compare <- compare %>% 
  rowwise() %>% 
  do(data.frame(., SplineDist(.$npoint, .$method))) %>% 
  ungroup()

ggplot(compare, aes(npoint, timeelas, colour = method)) +
  geom_line() +
  geom_point(aes(shape = method))

ggplot(compare, aes(npoint, perierror, colour = method)) +
  geom_line() +
  geom_point(aes(shape = method))

ggplot(compare, aes(npoint, itererror, colour = method)) +
  geom_line() +
  geom_point(aes(shape = method))
