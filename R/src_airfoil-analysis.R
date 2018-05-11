#============================#
# Airfoil Analysis
# Alwin Wang
#----------------------------#

#--- Generate single long_wall ----
AirfoilLongWall <- function(wallmsh) {
  # Determine trailing edge (most right point)
  te = wallmsh[which.max(wallmsh$x),1:2]
  # Determine leading edge point
  le <-  wallmsh %>% 
    mutate(dist = sqrt((x - te$x)^2 + (y - te$y)^2)) %>%
    arrange(-dist)
  le = le[1,1:2]
  # Determine centre point
  cp = (le + te)/2
  tetheta = atan2(te$y - cp$y, te$x - cp$x)
  # Number files
  wallmsh$wnum = 1:nrow(wallmsh)
  # Combine into long_wall
  long_wall <- wallmsh %>%
    mutate(theta = atan2(y - cp$y, x - cp$x) - tetheta) %>%
    mutate(theta = theta + ifelse(theta < 0, 2*pi, 0)) %>% # CANNOT use sign(theta) since theta can be zero
    arrange(theta, wnum)
  # Patch TE 
  # Note: It would be nice to decide if the first OR the second row should be used based on tail(...)
  long_wall$theta[1] = long_wall$theta[1] + 2*pi
  long_wall <- long_wall  %>%
    arrange(theta, wnum) %>%
    mutate(up = theta <= pi)
  # Patch LE if necessary
  lepatch <- filter(long_wall, x == le$x, y == le$y)
  if (nrow(lepatch) == 1) {
    # Add an extra LE row in 
    lepatch$bnum = NA; lepatch$wnum = NA
    lepatch$up = FALSE
    long_wall <- rbind(long_wall, lepatch) %>% arrange(theta)
  } else {
    long_wall$up <- ifelse(long_wall$wnum == max(lepatch$wnum), FALSE, long_wall$up)
  }
  # Add helpful variables
  long_wall <- mutate(long_wall, wall = !is.na(wnum))
  # Renumber long_wall
  long_wall$snum = 1:nrow(long_wall)
  # Return result
  return(long_wall)
}

#--- Airfoil Spline ----
# Calculations on closed simpled looped splines
AirfoilSpline <- function(long_wall, x = "x", y = "y", theta = "theta") {
  # Take only columns of interest, maybe consider making theta optional
  data = long_wall[,colnames(long_wall) %in% c(x, y, theta)]
  originalnrow = nrow(long_wall)
  # Remove duplicate rows or else cubic spline will give NaN
  # Note: spline is closed by repeated the first point, need a robust method for TE
  data = unique(data)
  # Rename columns of interest to be x and y
  oldnames <- colnames(data); names = oldnames
  names[names == x] = "x"
  names[names == y] = "y"
  names[names == theta] = "theta"
  colnames(data) <- names
  # Iterate to find the best spline distance
  data$s = c(0, EucDist(data$x, data$y)[-1])
  data$s = cumsum(data$s)
  # Initialise variable
  error = 1
  # Spline Length, s = int sqrt(dx/dt^2 + dy/dt^2) dt
  arclength <- function(tvec) {
    integral(function(t) sqrt(ppval(dcsx, t)^2 + ppval(dcsy, t)^2), tvec[1], tvec[2])}
  # Loop
  while (abs(error) > 0.001) {
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
    s <- pbapply(s, 1, arclength)
    # Determine cumulative spline distance
    s <- c(0, cumsum(s))
    # Report back error (loop variable)
    error = sum(data$s - s)
    print(error)
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
  # Combine data back with long_wall
  colnames(data) <- c(x, y, theta, colnames(data)[4:ncol(data)])
  long_wall <- full_join(long_wall, data, by = c(x, y, theta))
  # Check the number of rows has not increased
  if (originalnrow != nrow(long_wall)) {
    warning("Merge Failed")}
  # Return result
  return(long_wall)   # Data.frame
}

#--- Airfoil Offset ----
# Calculate distances from the surfaces (maybe requires analysis of one airfoil mesh)
AirfoilOffset <- function(long, totdist = 0.01, nsteps = 5, varh = FALSE) {
  # Use the cross product to determine the outward normal
  offset <- long$threaddata %>%
    filter(wall) %>%
    mutate(dirx = dyds*1 - 0*0,
           diry = -(dxds*1 - 0*0),
           dirdist = sqrt(dirx^2 + diry^2),
           dirx = dirx/dirdist,
           diry = diry/dirdist)
  # Repeat for 0:nsteps
  offset <- slice(offset, rep(1:n(), each = nsteps + 1))
  offset$nstep = rep(0:nsteps, length.out = nrow(offset)) # Length.out since nrow*5 already
  # Determine each of the distances
  offset <- offset %>%
    mutate(
      offseth = ifelse(is.na(aveh), totdist, ifelse(varh, aveh, totdist)),
      x = x + offseth*dirx*nstep/nsteps,
      y = y + offseth*diry*nstep/nsteps) %>%
    mutate(wall = ifelse(nstep == 0, TRUE, FALSE))
  # There are some duplicates in offset, I should remove them with UniLeftJoin!!
  # Plot
  # ggplot(offset, aes(x, y, colour = s, group = snum)) +
  #   geom_path() +
  #   geom_point(shape = 'o') +
  #   coord_fixed()
  # Return
  return(offset)   # Data.frame
}
