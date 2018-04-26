#============================
# Numerical Methods
# Alwin Wang
#----------------------------

#---- Splines ----
# Calculations on closed simpled looped splines
CalcSpline <- function(data) {
  t = atan2(data$y, data$x)
  t = t + ifelse(t < 0, 2*pi, 0)
  csx <- cubicspline(t, data$x)
  csy <- cubicspline(t, data$y)
  dcsx = csx; dcsx$coefs = t(apply(csx$coefs, 1, polyder))
  dcsy = csy; dcsy$coefs = t(apply(csy$coefs, 1, polyder))
  ds <- function(t) sqrt(ppval(dcsx, t)^2 + ppval(dcsy, t)^2)
  s = integral(ds, t[1], t[length(t)])
}