#============================
# Airfoil Analysis
# Alwin Wang
#----------------------------

#--- Generate long_bndry ----
# Create long_bndry
AirfoilLongBndry <- function(bndry, wallmsh) {
  # Determine trailing edge
  te = bndry[1,]
  # Determine leading edge point
  le = rbind(bndry, wallmsh) %>% 
    mutate(dist = sqrt((x - te$x)^2 + (y - te$y)^2)) %>%
    arrange(-dist)
  le = le[1,1:2]
  # Determine centre point
  cp = (le + te)/2
  tetheta = atan2(te$y - cp$y, te$x - cp$x)
  # Number files
  bndry$bnum = 1:nrow(bndry)
  bndry$wnum = NA
  wallmsh$bnum = NA
  wallmsh$wnum = 1:nrow(wallmsh)
  # Combine into long_bndry
  long_bndry <- rbind(bndry, wallmsh) %>%
    mutate(theta = atan2(y - cp$y, x - cp$x) - tetheta) %>%
    mutate(theta = theta + ifelse(theta < 0, 2*pi, 0)) # CANNOT use sign(theta) since theta can be zero
  # Assume that the boundary was closed at the TE. Thus, add 2pi to the "first" TE
  long_bndry$theta[1] = long_bndry$theta[1] + 2*pi
  long_bndry %<>% arrange(theta) %>%
    mutate(up = theta <= pi,
           snum = row_number())
  # Check that the trailing edge closes 
  # NOTE: indentical is too strong condition, all_equal could have been used
  if (!all_equal(long_bndry[1, 1:2], long_bndry[nrow(long_bndry), 1:2])) {
    warning("Trailing edge not closed")}
  # Patch LE if necessary
  lepatch <- filter(long_bndry, x == le$x, y == le$y)
  if (nrow(lepatch) == 1) {
    # Add an extra LE row in 
    lepatch$bnum = NA; lepatch$wnum = NA
    lepatch$up = FALSE
    long_bndry <- rbind(long_bndry, lepatch) %>% arrange(theta)
    # Renumber long_bndry
    long_bndry$snum = 1:nrow(long_bndry)
  } else {
    long_bndry$up <- ifelse(long_bndry$snum == max(lepatch$snum), FALSE, long_bndry$up)
  }
  # Add helpful variables
  long_bndry <- mutate(long_bndry, wall = !is.na(wnum), bndry = !is.na(bnum))
  # Plot
  ggplot(long_bndry, aes(x, y, colour = wall)) + #geom_path() + 
    geom_polygon(aes(group = wall, linetype = wall), fill = NA) + 
    # geom_point(aes(size = !is.na(wnum))) + 
    # scale_colour_gradientn(colours = myPalette(100)) +
    # coord_cartesian(xlim = c(0.599, 0.60478), ylim = c(-0.045, -0.04))
    coord_cartesian(xlim = c(0.60476, 0.60478), ylim = c(-0.042297, -0.042288))
  #--- RESULT ----
  # The boundary and wall mesh files should NOT be joined
  # This is because they do not coincide at the trailing edge
  # Return result
  return(long_bndry)
}



#--- Airfoil Spline ----
# Calculations on closed simpled looped splines
AirfoilSpline <- function(long_bndry, x = "x", y = "y") {
  # Take only columns of interest
  data = long_bndry[,colnames(long_bndry) %in% c(x, y)]
  # Remove duplicate rows or else cubic spline will give NaN
  data = unique(data)
  # Rename columns of interest to be x and y
  oldnames <- colnames(data); names = oldnames
  names[names == x] = "x"
  names[names == y] = "y"
  colnames(data) <- names
  # Iterate to find the best spline distance
  data$s = c(0, EucDist(data$x, data$y)[-1])
  data$s = cumsum(data$s)
  # Initialise variable
  error = 1
  # Spline Length, s = int sqrt(dx/dt^2 + dy/dt^2) dt
  length2 <- function(tvec) {
    integral(function(t) sqrt(ppval(dcsx, t)^2 + ppval(dcsy, t)^2), tvec[1], tvec[2])}
  # Loop
  while (abs(error) > 0.01) {
    # Cublic spline
    csx <- cubicspline(data$s, data$x)
    csy <- cubicspline(data$s, data$y)
    # Derivative of cubic splines
    dcsx = CubicSplineCalc(csx, -1)
    dcsy = CubicSplineCalc(csy, -1)
    # Determine the intervals for arc lengths (quicker than cumulative lengths)
    s <- cbind(data$s, lead(data$s))[-nrow(data),]
    # Determine interval spline distances
    cat("Determining spline distance\n")
    s <- pbapply(s, 1, length2)
    # Determine cumulative spline distance
    s <- c(0, cumsum(s))
    # Report back error (loop variable)
    error = sum(data$s - s)
    # Update the loop
    data$s = s
  }
  # Cublic spline
  csx <- cubicspline(data$s, data$x)
  csy <- cubicspline(data$s, data$y)
  # Derivative of cubic splines
  dcsx = CubicSplineCalc(csx, -1)
  dcsy = CubicSplineCalc(csy, -1)
  data$dxds <- ppval(dcsx, data$s)
  data$dyds <- ppval(dcsy, data$s)
  data$dydx <- data$dyds/data$dxds
  # Combine data back with long_bndry
  colnames(data) <- c(x, y, colnames(data)[3:ncol(data)])
  long_bndry <- full_join(long_bndry, data, by = c(x, y))
  return(long_bndry)
}