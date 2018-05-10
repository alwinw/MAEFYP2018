#============================
# Numerical Methods
# Alwin Wang
#----------------------------

#--- Required Functions ----
# Heavyside function (step function)
heav <- function(t) ifelse(t>0,1,0)

# Distance function
EucDist <- function(x, y) sqrt((x - lag(x))^2 + (y - lag(y))^2)

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

#--- Interpolation Function ----
interpolate <- function(omesh, imesh, 
                        onames = c("x", "y", "z"), inames = c("x", "y"),
                        wallsplit = TRUE) {
  # Selet data based on onames for omesh
  
  # Check for NAs
  
  # Split datat into wall and non-wall
  
  # Interpolate onto imesh using duplicate = "strip"
  
  # Recombine wall and non wall if wallsplit == TRUE
  
  # Check the interpolation accuracy by interpolating back onto omesh and save column % error
  
  # Add original column names
  
}