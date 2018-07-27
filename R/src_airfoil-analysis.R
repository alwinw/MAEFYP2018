#============================#
# Airfoil Analysis
# Alwin Wang
#----------------------------#

#--- Generate single long_wall ----
AirfoilLongWall <- function(wallmsh) {
  # Determine trailing edge (most right point)
  te = wallmsh[which.max(wallmsh$x),1:2]
  # Determine leading edge point (furthest point from the TE)
  le <-  wallmsh %>% 
    mutate(dist = sqrt((x - te$x)^2 + (y - te$y)^2)) %>%
    arrange(-dist)
  le = le[1,1:2]
  # Determine centre point and angle to the te
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
    # mutate(theta = 2*pi - theta) %>% # start from TE -> lower surface -> LE -> upper surface -> TE
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
  data$dydxlen <- sqrt(data$dxds^2 + data$dyds^2)
  # Combine data back with long_wall
  colnames(data) <- c(x, y, theta, colnames(data)[4:ncol(data)])
  long_wall <- full_join(long_wall, data, by = c(x, y, theta))
  # Check the number of rows has not increased
  if (originalnrow != nrow(long_wall)) {
    warning("Merge Failed")}
  # Return result
  return(long_wall)   # Data.frame
}

AirfoilSplineCheck <- function(long_wall) {
  splinecheck <- long_wall %>%
    select(x, y, nxG, nyG, up, dxds, dydx, dydx) %>%
    mutate(dydxG = -nxG/nyG) %>%
    mutate(error = abs((dydxG-dydx)/dydx))
  if (cor(splinecheck$dydx, splinecheck$dydxG) < 0.99)
    warning("Low correlation between spline and wallgrad")
  if (max(splinecheck$error) > 0.01)
    warning(paste(sum(splinecheck$error > 0.01), "spline abs error > 1%"))
  if (mean(splinecheck$error) > 0.0015)
    warning("Mean spline abs error > 0.15%")
  # ave(gradient1, gradient2) =/= ave(Dely1, Dely2)/ave(Delx1, Delx2)
  splinecheck <- splinecheck %>%
    group_by(x, y, up) %>%
    mutate(nxG = mean(nxG), nyG = mean(nyG)) %>%
    ungroup() %>%
    mutate(dydxGm = -nxG/nyG) %>%
    mutate(errorm = abs((dydxGm - dydx)/dydx))
  # cor(splinecheck$dydx, splinecheck$dydxG) <= cor(splinecheck$dydx, splinecheck$dydxGm)
  # max(splinecheck$error)  >= max(splinecheck$errorm)
  # mean(splinecheck$error) >= mean(splinecheck$errorm)
  
  long_wall <- long_wall %>%
    group_by(x, y, up) %>%
    mutate(nxG = mean(nxG), nyG = mean(nyG), areaG = mean(areaG)) %>%
    ungroup()
  
  return(long_wall)
}

#--- Airfoil Offset ----
# Calculate distances from the surfaces (maybe requires analysis of one airfoil mesh)
AirfoilOffset <- function(long, totdist = 0.008, nsteps = 5, varh = TRUE, scale) {
  # Use the cross product to determine the outward normal
  offset <- long$threaddata %>%
    filter(wall) %>%
    group_by(s) %>% # mutate(temp_uni = ifelse(is.na(aveh), totdist, aveh))
    top_n(1, aveh) %>%
    ungroup() %>%
    mutate(dirx = dyds*1 - 0*0,
           diry = -(dxds*1 - 0*0),
           dirdist = sqrt(dirx^2 + diry^2),
           dirx = dirx/dirdist,
           diry = diry/dirdist)
  # Add a number to offset
  offset$onum <- 1:nrow(offset)
  # Repeat for 0:nsteps
  offset <- slice(offset, rep(1:n(), each = nsteps + 1))
  offset$nstep = rep(0:nsteps, length.out = nrow(offset)) # Length.out since nrow*5 already
  # Determine each of the distances
  offset <- offset %>%
    mutate(
      offseth = ifelse(is.na(aveh) | !varh, totdist, aveh)*scale,
      norm = offseth*nstep/nsteps,
      x = x + dirx*norm,
      y = y + diry*norm) %>%
    mutate(wall = ifelse(nstep == 0, TRUE, FALSE),
           wnum = ifelse(wall, wnum, NA))
  # There are some duplicates in offset, I should remove them with UniLeftJoin!!
  # Plot
  # ggplot(offset, aes(x, y, colour = s, group = snum)) +
  #   geom_path() +
  #   geom_point(shape = 'o') +
  #   coord_fixed()
  # Return
  return(offset)   # Data.frame
}

# Calculate actual enum of offset points
AirfoilOffsetEnum <- function(long, localnum = 2, returnList = FALSE) {
  long$offset$enum_ori = long$offset$enum
  # Data
  poly_df <- long$threaddata %>% 
    filter(local <= localnum, seshnode) %>%
    arrange(enum, ncorner) %>%
    select(x, y, enum)
  pts_df <- long$offset %>%
    select(x, y)
  # Plot
  # ggplot() + 
  #   geom_polygon(aes(x, y, group = enum), 
  #     poly_df, fill = NA, colour = "black") +
  #   geom_point(aes(x, y), pts_df, shape = 'O') + 
  #   coord_fixed()
  # Split the dataframe into a list based on enum and then remove enum from df in the list
  poly_list <- split(poly_df, poly_df$enum)
  # Convert the list to Polygon, then create a Polygons object
  poly_sp <- sapply(poly_list, function(poly){
    Polygons(list(Polygon(poly[, c("x", "y")])), ID = poly[1, "enum"])
  })
  # polygonsp <- Polygons(polygonsp, ID = 1)
  poly_sp <- SpatialPolygons(poly_sp)
  # plot(poly_sp)
  # plot(poly_sp, col=poly_sp@plotOrder)
  # Convert points to coordinates
  pts_ps <- pts_df
  coordinates(pts_ps) <- ~x+y
  # points(pts_ps$x, pts_ps$y)
  # Determine polygons points are in
  pts_return <- over(pts_ps, poly_sp, returnList = returnList)
  # points(pts_ps$x, pts_ps$y, col = pts_return)
  if (returnList) {
    pts_return <- lapply(pts_return, function(pt) {
      if(length(pt) == 0) data.frame(enum = NA, ptin = NA)
      else data.frame(enum = pt, ptin = length(pt))
    })
    # Recombine
    pts_list <- split(pts_df, rownames(pts_df))
    # DO SOMETHING FANCY HERE
    temp <- mapply(c, pts_list, pts_return, SIMPLIFY = FALSE)
    temp <- lapply(temp, function(pt) {
      data.frame(x = pt$x, y = pt$y, enum = pt$enum, ptin = pt$ptin)
    })
    temp <- bind_rows(temp)
    temp$enum <- unique(poly_df$enum)[temp$enum]
    warning("Not finished yet")
  } else {
    # output <- cbind(pts_df, data.frame(enum = pts_sp))
    long$offset$enum = unique(poly_df$enum)[pts_return]
  }
  # Return
  return(long)
}

#--- Airfoil Coordinate Transform ----
AirfoilTransform <- function(long, localnum = 2, extrap = 0.05)  {
  # Limits of splines
  lim_airfoil <- data.frame(
    te_up = min(long$walldata$s),
    le = long$walldata[long$walldata$theta == pi,]$s[1],
    te_lo = max(long$walldata$s)) %>%
    mutate(min = te_up - extrap, max = te_lo + extrap)
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
  # Find the minimum distance for all points
  pts_airfoil_upper <- long$threaddata %>% 
    filter(local <= localnum) %>%
    select(x, y, s) %>%
    unique() %>% 
    rowwise() %>%
    do(data.frame(., MinS(.$x, .$y, csx, csy, dcsx, dcsy,
                          lim_airfoil$te_up, lim_airfoil$te_lo)))
  # Plots
  # ggplot(pts, aes(s)) + geom_density()
  # ggplot(pts, aes(prec)) + geom_density()
  # ggplot(pts, aes(dotprod)) + geom_density()
  # Separate out the wake nodes into upper and lower
  pts <- pts_airfoil_upper %>%
    mutate(coord = ifelse(abs(dotprod) < 1e-5, "Airfoil", NA))
  pts_upper <- pts %>%
    filter(is.na(coord)) %>%
    select(x, y, s) %>%
    do(data.frame(., MinS(.$x, .$y, csx, csy, dcsx, dcsy,
                          lim_airfoil$min, lim_airfoil$te_up))) %>%
    mutate(coord = ifelse(abs(dotprod) < 1e-5 & crossprod > 0, 
                          "Upper", NA))
  pts_lower <- pts %>%
    filter(is.na(coord)) %>%
    select(x, y, s) %>%
    do(data.frame(., MinS(.$x, .$y, csx, csy, dcsx, dcsy,
                          lim_airfoil$te_lo, lim_airfoil$max))) %>%
    mutate(coord = ifelse(abs(dotprod) < 1e-5 & crossprod > 0,
                          "Lower", NA))
  # Clean up
  pts <- rbind(pts, pts_upper, pts_lower) %>%
    filter(!is.na(coord),
           stream != lim_airfoil$min,
           stream != lim_airfoil$max)
  # Plots
  long$airfoilextrap <- data.frame(
   s = seq(lim_airfoil$te_up - extrap, lim_airfoil$te_lo + extrap, length.out = 100)) %>%
   mutate(x = ppval(csx, s), y = ppval(csy, s))
  # ggplot() +
  #   geom_path(aes(x, y), long$airfoilextrap) +
  #   geom_point(aes(x, y, colour = coord), pts) +
  #   geom_path(aes(x, y, group = snum), long$offset) +
  #   coord_fixed(xlim = c(0.4, 0.75))
  #  # coord_fixed()
  
  # Join with original data
  long$threaddata <-  LongJoin(long$threaddata, select(pts, -prec, dotprod, crossprod, dist),
           unicols = c("x", "y", "s", "coord"))
  # Return the output
  return(long)
}


