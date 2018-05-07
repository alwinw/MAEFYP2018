#============================
# Numerical Methods
# Alwin Wang
#----------------------------

#--- Required Functions ----
# Heavyside function (step function)
heav <- function(t) ifelse(t>0,1,0)

#--- Cubic Spline Calculus ----
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
  # Plots
  # plot <- data.frame(x = seq(cs$breaks[1], cs$breaks[length(cs$breaks)], length = 100))
  # plot$y = ppval(cs, plot$x)
  # plot$cy = ppval(ccs, plot$x)
  # ggplot(plot, aes(x = x)) +
  #   geom_path(aes(y = y)) + geom_point(aes(y = y)) +
  #   geom_path(aes(y = cy)) + geom_point(aes(y = cy))
  return(ccs)
}

#--- Splines ----
# Calculations on closed simpled looped splines
# e.g. rawdata = long_bndry; colnames(rawdata) <- c("hi", "bye", "up", "wnum", "wsnum", "snum"); x = "hi"; y = "bye"
CalcSpline <- function(rawdata, x = "x", y = "y") {
  # Take only columns of interest
  data = rawdata[,colnames(rawdata) %in% c(x, y)]
  # Remove duplicate rows or else cubic spline will give NaN
  data = unique(data)
  # Rename columns of interest to be x and y
  oldnames <- colnames(data); names = oldnames
  names[names == x] = "x"
  names[names == y] = "y"
  colnames(data) <- names
  # Use the path distance from the first point as the parametric variable, t
  data$t = c(0, sqrt((data$x - lag(data$x))^2 + (data$y - lag(data$y))^2)[-1])
  data$t = cumsum(data$t)
  # Cublic spline
  csx <- cubicspline(data$t, data$x)
  csy <- cubicspline(data$t, data$y)
  # Derivative of cubic splines
  dcsx = CubicSplineCalc(csx, -1)
  dcsy = CubicSplineCalc(csy, -1)
  # Integral of cubic splines
  # icsx = CubicSplineCalc(csx, 1)
  # icsy = CubicSplineCalc(csy, 1)
  # Spline Length, s = int sqrt(dx/dt^2 + dy/dt^2) dt
  length2 <- function(tvec) integral(function(t) sqrt(ppval(dcsx, t)^2 + ppval(dcsy, t)^2), tvec[1], tvec[2])
  s <- cbind(data$t, lead(data$t))[-nrow(data),]
  s <- apply(s, 1, length2)
  s <- c(0, cumsum(s))
  # Derivative dy/dx = dy/dt * dt/dx
  # !! Maybe consider refitting the cubic polynomials on s rather than t
  dxdt <- ppval(dcsx, data$t)
  dydt <- ppval(dcsy, data$t)
  dydx <- dydt/dxdt
  # Plot
  # plot <- data.frame(t = data$t, dydx = dydx, x = data$x, y = data$y)
  # myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  # ggplot(plot, aes(x = x, y = y, colour = ifelse(abs(dydx)>0.5,0.5*sign(dydx),dydx))) + geom_path() +
  #   scale_colour_gradientn(colours = myPalette(100)) + geom_point(shape=as.logical(sign(dydx)/2+1/2)) 
  # Perhaps add a basic check of dy/dx vs y2-y1/x2-x1 and s vs t
  # Combine results with data
  data$s = s
  data$dydx = dydx
  data$t <- NULL
  # Combine data back with rawdata
  colnames(data) <- c(x, y, colnames(data)[3:ncol(data)])
  rawdata <- full_join(rawdata, data, by = c(x, y))
  return(rawdata)
}
