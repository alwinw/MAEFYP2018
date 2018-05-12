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
      norm = offseth*nstep/nsteps,
      x = x + dirx*t,
      y = y + diry*t) %>%
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

#--- Airfoil Coordinate Transform ----
AirfoilTransform <- function(long, extrap = 0.2)  {
  # Limits of splines
  lim_airfoil <- data.frame(
    te_up = min(long$walldata$s),
    le = long$walldata[long$walldata$theta == pi,]$s[1],
    te_lo = max(long$walldata$s))
  # Determine unique wall points for cubic spline
  uni_wall <- long$walldata %>%
    select(x, y, s) %>%
    unique()
  # Determine spline polynomials
  csx <- cubicspline(uni_wall$s, uni_wall$x)
  csy <- cubicspline(uni_wall$s, uni_wall$y)
  # Derivative of cubic splines
  dcsx = CubicSplineCalc(csx, -1)
  dcsy = CubicSplineCalc(csy, -1)
  # Function to find minimum distance
  MinS <- function(x, y) {
    # Determine minimum distance (bounded!)
    min.out <- fminbnd(
      function(s) {sqrt((x - ppval(csx, s))^2 + (y - ppval(csy, s))^2)},
      lim_airfoil$te_up, lim_airfoil$te_lo)
    # Determine dot product (if time, also normalise it)
    dotprod = 
      (x - ppval(csx, min.out$xmin))*ppval(dcsx, min.out$xmin) +
      (y - ppval(csy, min.out$xmin))*ppval(dcsy, min.out$xmin)
    # Create and return output
    return(data.frame(
      s.min = min.out$xmin,
      norm.min = min.out$fmin,
      prec.min = min.out$estim.prec,
      dotprod.min = dotprod))
  }
  # Find the minimum distance for all points
  pts <- long$threaddata %>% 
    filter(local <= 3) %>%
    # select(x, y, s) %>% 
    unique() %>% 
    rowwise() %>%
    do(data.frame(., MinS(.$x, .$y)))
  # Plots
  # ggplot(pts, aes(s)) + geom_density()
  # ggplot(pts, aes(prec.min)) + geom_density()
  # ggplot(pts, aes(dotprod.min)) + geom_density()
  ggplot() +
    geom_point(aes(x, y, colour = abs(dotprod.min) > 1e-5), pts) +
    geom_path(aes(x, y, group = snum), long$offset) +
    coord_fixed(xlim = c(0.4, 0.75))
  ggplot(pts %>% filter(abs(dotprod.min) < 1e-5, !is.na(nnum)) %>% arrange(enum, ncorner), 
         aes(s.min, norm.min, group = enum, colour = enum)) +
    geom_polygon(fill = NA) +
    geom_point(data = pts %>% filter(abs(dotprod.min) < 1e-5))
}



#--- Airfoil Coordinate Transform ----
AirfoilTransform <- function(long, extrap = 0.2) {
  
  # Determine unique wall points for cubic spline
  uni_wall <- long$walldata %>%
    select(x, y, s) %>%
    unique()
  # Determine spline polynomials
  csx <- cubicspline(uni_wall$s, uni_wall$x)
  csy <- cubicspline(uni_wall$s, uni_wall$y)
  # Derivative of cubic splines
  dcsx = CubicSplineCalc(csx, -1)
  dcsy = CubicSplineCalc(csy, -1)
  # Points to coordinate transform (not necessarily threaddata!!)
  long_ori <- long$threaddata %>%
    filter(!is.na(nnum), local <= 2) %>%
    select(x, y) %>%
    unique() %>%
    mutate(onum = row_number())
  # Limits
  lim <- data.frame(lim_s = c(lim_airfoil$te_up, lim_airfoil$le))
  lim$id = c("lo", "up")
  lim$lim_x = ppval(csx, lim$lim_s) 
  lim$lim_y = ppval(csy, lim$lim_s)
  lim$lim_dx = ppval(dcsx, lim$lim_s)
  lim$lim_dy = ppval(dcsy, lim$lim_s)
  # Determine if (x,y) lie in this coordiante transform
  long_lim <- slice(long_ori, rep(1:n(), each = nrow(lim)))
  long_lim <- cbind(long_lim, slice(lim, rep(1:n(), length.out = nrow(long_lim))))
  
  long_lim <- long_lim %>% 
    mutate(
      vecx = x - lim_x,
      vecy = y - lim_y,
      dotprod = vecx*lim_dx + vecy*lim_dy,
      crossprod = vecx*lim_dy - vecy*lim_dx,
      lobound = dotprod >= 0 & crossprod >= 0,
      upbound = dotprod <=0 & crossprod >= 0) %>%
    group_by(onum) %>%
    mutate(
      bounded = sum(lobound, upbound)) 
  # %>%
    # filter(id == "up")
  
  # ggplot(long_lim, aes(x, y, colour = bounded == 2)) + geom_point() + coord_fixed()
  
  lim <- data.frame(lim_s = seq(0, 2, length.out = 100))
  lim$lim_x = ppval(csx, lim$lim_s) 
  lim$lim_y = ppval(csy, lim$lim_s)
  lim$lim_dx = ppval(dcsx, lim$lim_s)
  lim$lim_dy = ppval(dcsy, lim$lim_s)
  
  long_lim <- slice(data.frame(x = 0, y = -1), rep(1:n(), each = nrow(lim)))
  long_lim <- cbind(long_lim, slice(lim, rep(1:n(), length.out = nrow(long_lim))))
  
  long_lim <- long_lim %>% 
    mutate(
      vecx = x - lim_x,
      vecy = y - lim_y,
      dotprod = vecx*lim_dx + vecy*lim_dy,
      crossprod = vecx*lim_dy - vecy*lim_dx,
      lobound = dotprod >= 0 & crossprod >= 0,
      upbound = dotprod <=0 & crossprod >= 0) %>%
    # group_by(onum) %>%
    mutate(
      bounded = sum(lobound, upbound)) 
  
  ggplot(long_lim, aes(lim_s, dotprod, colour = crossprod)) + 
    geom_point() +
    geom_path()
  
  ggplot(long_lim, aes(lim_s, dotprod^sign(crossprod), colour = crossprod)) + 
    geom_point() +
    geom_path()
}



