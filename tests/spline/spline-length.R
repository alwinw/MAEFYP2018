# Libraries
library("bezier")
library("pracma")
library("dplyr")
library("ggplot2")

# Generate random data
theta = seq(0, 2*pi, length.out = 50) + runif(50, -0.1, 0.1)
theta = sort(theta)
theta =  c(0, theta[theta >=0 & theta <= 2*pi], 2*pi)
data = data.frame(x = cos(theta), y = sin(theta))

ggplot(data, aes(x, y)) + geom_path() + geom_point()

# Euclidean Distance
EucDist <- function(data) {
  dist <- sqrt((data$x - lag(data$x))^2 + (data$y - lag(data$y))^2)
  dist <- sum(dist[2:length(dist)])
  return(dist)
}

# Cubic Spline Full Integration
CSFull <- function(data) {
  t = atan2(data$y, data$x)
  t = t + ifelse(t < 0, 2*pi, 0)
  csx <- cubicspline(t, data$x)
  csy <- cubicspline(t, data$y)
  dcsx = csx; dcsx$coefs = t(apply(csx$coefs, 1, polyder))
  dcsy = csy; dcsy$coefs = t(apply(csy$coefs, 1, polyder))
  ds <- function(t) sqrt(ppval(dcsx, t)^2 + ppval(dcsy, t)^2)
  s = integral(ds, t[1], t[length(t)])
}

# Do an iterative piecewise one, check each ?integral option

