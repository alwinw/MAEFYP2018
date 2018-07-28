#============================#
# Helper Functions
# Alwin Wang
#----------------------------#

#--- Airfoil Calculation                                          ----
#--- * Boundary Data                                              ----
# Load the boundary file (i.e. airfoil shape)
  # out: bndry = data.frame(x, y)
LoadBndry <- function(folder,
                      profile = "bndry_prf", extension = ".dat") {
  file            <- paste0(folder, profile, extension)
  bndry           <- read.table(file)
  colnames(bndry) <- c("x", "y")
  # out: bndry = data.frame(x, y)
  return(bndry)
}
#--- * Wall Mesh Data                                             ----
# Load wall grad file
# out: wallmesh = data.frame(x, y, nxG, nyG, areaG)
LoadWallGrad <- function(seshpath,
                         extension = ".wallgrad") {
  file               <- paste0(seshpath, extension)
  wallmesh           <- read.table(file, skip = 1)
  colnames(wallmesh) <- c("x", "y", "nxG", "nyG", "areaG")
  # out: wallmesh = data.frame(x, y, nxG, nyG, areaG)
  return(wallmesh)
}
# Determine the LE, TE, spline length, etc of the airfoil
# out:long_wall = data.frame(<wallmesh>, wnum, theta, up, wall, snum)
AirfoilLongWall <- function(wallmesh, 
                            dir = "clockwise") {
  # Determine LE and TE edges
  te = wallmesh[which.max(wallmesh$x), 1:2]                       # Most far right point
  le <- wallmesh %>%                                              # Furthest point from te
    mutate(dist = sqrt((x-te$x)^2 + (y-te$y)^2)) %>% arrange(-dist)
  le = le[1, 1:2]
  # Determine theta and sort
  cp = (le+te)/2
  tetheta = atan2(te$y-cp$y, te$x-cp$x)
  wallmesh$wnum = 1:nrow(wallmesh)
  long_wall <- wallmesh %>%
    mutate(theta = atan2(y - cp$y, x - cp$x) - tetheta) %>%
    mutate(theta = theta + ifelse(theta < 0, 2*pi, 0)) %>%        # CANNOT use sign(theta) since theta can be zero
    arrange(theta, wnum)
  # Patch TE 
  long_wall$theta[1] = long_wall$theta[1] + 2*pi
  # Order by theta
  if (dir == "clockwise") {
    long_wall <- long_wall  %>%
      mutate(theta = 2*pi - theta) %>%                              # TE -> lower -> LE -> upper -> TE
      arrange(theta, wnum) %>%
      mutate(up = theta <= pi)
  } else {
    long_wall <- long_wall  %>%
      arrange(theta, wnum) %>%
      mutate(up = theta <= pi)
  }
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
  long_wall <- mutate(long_wall, wall = !is.na(wnum))             # Add helpful variables
  long_wall$snum = 1:nrow(long_wall)
  # out:long_wall = data.frame(<wallmesh>, wnum, theta, up, wall, snum)
  return(long_wall)
}
# Determine airfoil spline length
# out: long_wall = data.frame(<long_wall>, s, dxds, dyds, dydx, dydxlen)
AirfoilSpline <- function(long_wall, 
                          x = "x", y = "y", theta = "theta") {
  # Take only columns of interest, maybe consider making theta optional
  data = long_wall[,colnames(long_wall) %in% c(x, y, theta)]
  originalnrow = nrow(long_wall)
  # Note: spline is closed by repeated the first point, need a robust method for TE
  data = unique(data)                                             # Remove duplicate rows
  oldnames <- colnames(data); names = oldnames                    # Rename columns of interest to be x and y
  names[names == x] = "x"
  names[names == y] = "y"
  names[names == theta] = "theta"
  colnames(data) <- names
  # Iterate to find the best spline distance
  data$s = c(0, LagDist(data$x, data$y)[-1])
  data$s = cumsum(data$s)
  error = 1                                                       # Initialise variable
  while (abs(error) > 0.001) {                                    # Loop for required error
    csx <- cubicspline(data$s, data$x)                            # Cublic spline
    csy <- cubicspline(data$s, data$y)
    dcsx = CubicSplineCalc(csx, -1)                               # Derivative of cubic splines
    dcsy = CubicSplineCalc(csy, -1)
    cat("Determining spline distance\n")
    s <- cbind(data$s, lead(data$s))[-nrow(data),]                # Interval (s1, s2)
    s <- pbapply(s, 1, CubicSplineArcLength, dcsx, dcsy)          # Arc length of interval
    s <- c(0, cumsum(s))                                          # Cumulative arc lenth
    error  = sum(data$s - s)                                      # Error
    data$s = s                                                    # Update loop
    print(error)
  }
  csx <- cubicspline(data$s, data$x)                              # Cublic spline
  csy <- cubicspline(data$s, data$y)
  dcsx = CubicSplineCalc(csx, -1)                                 # Derivative of cubic splines
  dcsy = CubicSplineCalc(csy, -1)
  data$dxds    <- ppval(dcsx, data$s)                             # Use cublic spline to determine derivatives
  data$dyds <- ppval(dcsy, data$s)
  data$dydx <- data$dyds/data$dxds
  data$dydxlen <- sqrt(data$dxds^2 + data$dyds^2)
  # Combine data back with long_wall
  colnames(data) <- c(x, y, theta, colnames(data)[4:ncol(data)])
  long_wall <- full_join(long_wall, data, by = c(x, y, theta))
  if (originalnrow != nrow(long_wall)) {                          # Check the number of rows is correct
    warning("Merge Failed")}
  # out: long_wall = data.frame(<long_wall>, s, dxds, dyds, dydx, dydxlen)
  return(long_wall)
}

#--- Session and Mesh Calculation                                 ----
#--- * Session Data                                               ----
# Load a single keyword from list of keywords
LoadKeyword <- function(keyword, 
                        filelines, file) {
  num <- grep(keyword[1], filelines)                              # Find the line numbers the keyword appears in
  # keyword... keyword
  table <- read.table(
    file = file,
    skip = num[1],
    nrow = num[2]-num[1]-1,
    stringsAsFactors = FALSE)
  colnames(table) <- keyword[2:length(keyword)]                   # Set talble colname names
  table$junk <- NULL                                              # Remove any junk defined in "keyword"
  return(table)
}
# Load a list of keywords and colnames
# out: session = list(nodes, elements, surfaces, curves)
LoadSeshFileKeywords <- function(seshpath, 
                                 list_keyword = NULL) {
  # Determine keywords
  if (is.null(list_keyword)) {
    list_keyword <- list(                                         # Keywords in session to read
      c("NODES", "nnum", "x", "y", "z"),
      c("ELEMENTS", "enum", "shapetag", "n1", "n2", "n3", "n4", "junk"),
      c("SURFACES", "snum", "element", "side", "bctag", "bc", "junk"),
      c("CURVES", "cnum", "element", "side", "curvetag", "curvedata", "junk"))
  }
  # Read keywords from file
  file      <-  paste0(seshpath,".sesh")                          # Session file name
  filelines <- readLines(file)                                    # Read the session file to grep lines later
  session   <- lapply(list_keyword, LoadKeyword, filelines, file) # For each keyword load the associated data
  names(session) <- sapply(list_keyword,                          # Set the names for the list
                           function(key) tolower(key[1]))
  # out: session = list(nodes, elements, surfaces, curves)
  return(session)
}
# Determine the mesh nodes from the session data
# out: long_sesh = tibble(enum, ncorner, nnum, x, y, z, elabx, elaby, area
LongSesh <-function(session) {
  # Check that all elements are quadrangles
  if (as.logical(sum(session$elements$shapetag != "<Q>"))) {
    warning("Not all elements are Quadrangles")
  }
  # Manipulate data
  long_sesh <- session$elements %>%                           # data.frame(enum, shapetag, n1, n2, n3, n4)
    select(-shapetag) %>%                                         # Remove shape tag since checked shape already
    gather(ncorner, nnum, -enum) %>%                              # Gather nodes 1-4 into a single column
    left_join(., session$nodes, "nnum") %>%                       # Combine with nodes to get (x,y) for each element corner
    group_by(enum) %>%                                            # Order data
    mutate(                                                       # Add columns for labels
      elabx = mean(x),
      elaby = mean(y)) %>%
    rbind(.,                                                      # Close the quadrilateratl (without effecting mean) 
          mutate(filter(., ncorner=="n1"), ncorner = "n5")) %>%
    arrange(enum, ncorner) %>% 
    mutate(area = x*lead(y) - lead(x)*y) %>%                      # Calculate area
    mutate(area = 1/2*abs(sum(ifelse(is.na(area), 0, area)))) %>%
    filter(ncorner != "n5") %>%                                   # Remove extra ncorner
    ungroup()
  # out: long_sesh = tibble(enum, ncorner, nnum, x, y, z, elabx, elaby, area)
  return(long_sesh)
}
#--- * Mesh Data                                                  ----
# Load mesh file
# out: mesh = data.frame(x, y, enum, jnum)
LoadMesh <- function(seshpath,
                     extension = ".mshi") {
  file           <- paste0(seshpath, extension)
  mesh           <- read.table(file, skip = 1)
  colnames(mesh) <- c("x", "y", "enum", "jnum")
  # out: mesh = data.frame(x, y, enum, jnum)
  return(mesh)
}
# Convert to long mesh
# out: long_mesh = data.frame(x, y, enum, jnum, ncorner, nnum, z, elabx, elaby, area, node) 
LongMesh <- function(long_mesh, long_sesh) {
  long_mesh$mnum =  1:nrow(long_mesh)                              # Mesh number
  long_mesh     <- LongJoin(long_mesh, long_sesh)                  # Join with node numbers
  long_mesh$z   <- NULL                                            # Not useful column
  long_mesh$node = !is.na(long_mesh$nnum)                          # Handy variable
  # out: long_mesh = data.frame(x, y, enum, jnum, ncorner, nnum, z, elabx, elaby, area, node) 
  return(long_mesh)
}
#--- * Wall Data                                                  ----
LongWall <- function(long_wall, long_mesh) {
  # Non-node points only
  join_mesh <- long_mesh %>%
    filter(!node) %>%
    select(x, y, enum)
  join_wall <- left_join(long_wall, join_mesh, by = c("x", "y")) %>%
    arrange(wnum)
  # Path up NA enum (assume wnum was some logical order)
  join_wall$node <- is.na(join_wall$enum)
  join_wall$enum <- ifelse(is.na(join_wall$enum), lag (join_wall$enum), join_wall$enum)
  join_wall$enum <- ifelse(is.na(join_wall$enum), lead(join_wall$enum), join_wall$enum)
  # Double check
  if (length(unique(count(group_by(join_wall, enum))[,2])) != 1)
    warning("Unequal number of points per element")
  # Determine average element height
  area <- long_mesh %>%
    filter(node, enum %in% unique(join_wall$enum)) %>%
    select(enum, area) %>%
    arrange(enum) %>%
    unique(.)
  base <- join_wall %>%
    filter(node) %>%
    group_by(enum) %>%
    mutate(base = LagDist(x, y)) %>%
    filter(!is.na(base)) %>%
    select(enum, base) %>%
    arrange(enum) %>%
    unique(.)
  if (sum(abs(base$enum - area$enum)) > 0)
    warning("Unexplained error between area and base")
  elem <- cbind(area, base = base$base) %>%
    mutate(aveh = area/base,
           ar   = base/aveh)
  long_wall <- left_join(join_wall, elem, by = "enum") %>%
    arrange(snum)
  return(long_wall)
}


#--- Numerical Methods ----
# Heavyside function (step function)
heav <- function(t) ifelse(t>0,1,0)

# Distance function
LagDist <- function(x, y) sqrt((x - lag(x))^2 + (y - lag(y))^2)

# Function to find minimum distance
MinS <- function(x, y, csx, csy, dcsx, dcsy, lim_low, lim_up) {
  # Determine minimum distance (bounded!)
  min.out <- fminbnd(
    function(s) {sqrt((x - ppval(csx, s))^2 + (y - ppval(csy, s))^2)},
    lim_low, lim_up)
  # Determine vectors
  vecpt = data.frame(
    x = x - ppval(csx, min.out$xmin),
    y = y - ppval(csy, min.out$xmin))
  vecsf = data.frame(
    x = ppval(dcsx, min.out$xmin),
    y = ppval(dcsy, min.out$xmin))
  # Determin dot products
  dotprod = vecpt$x*vecsf$x + vecpt$y*vecsf$y
  crossprod = vecpt$x*vecsf$y - vecsf$x*vecpt$y
  dist = sqrt(vecpt$x^2 + vecpt$y^2) * sqrt(vecsf$x^2 + vecsf$y^2)
  # Note: if I do dotprod/dist, at the surface ~ 0/0 != 0 because of machine error
  # Create and return output
  return(data.frame(
    stream = min.out$xmin,
    norm = min.out$fmin,
    prec = min.out$estim.prec,
    dotprod = dotprod,
    crossprod = crossprod,
    dist = dist))
}

# Unique values
UniLeftJoin <- function(longdata, unicols = c("x", "y")) {
  longdata[!duplicated(longdata[,unicols]),]
}

# Custom left join
LongJoin <- function(left_data, right_data, unicols = NULL, wall = FALSE) {
  # Prepare data for joining
  join_cols <-                                                    # Determine columns to join over
    colnames(left_data)[colnames(left_data) %in% colnames(right_data)]
  if (is.null(unicols)) unicols = join_cols                       # Cols to determine duplicate rows
  uni_right_data <- UniLeftJoin(right_data, unicols)              # Remove duplicate rows
  # Join left_data and right_data
  join_data <- left_join(                                         # Join ALL columns and rows
    left_data, uni_right_data, by = join_cols) 
  if (wall) join_data <- mutate(join_data, wall = !is.na(wnum))   # Update wall column if wanted
  # Check for duplicates
  if(nrow(join_data) != nrow(left_data)) {                        # Check the number of rows has not increased
    warning(paste(
      "join_data", "has more rows than", 
      deparse(substitute(left_data)),
      "by", nrow(join_data) - nrow(left_data)))}
  # Check for missing values
  join_anti <- anti_join(                                         # Find rows not merged in
    uni_right_data, left_data, by = join_cols)
  if(nrow(join_anti) != 0) {
    warning(paste(
      deparse(substitute(right_data)), 
      "missing rows in", "join_data",
      "by", nrow(join_anti)))}
  # Return Result
  return(join_data)
}

#--- Cubic Spline Calculus ----
# Spline Length
CubicSplineArcLength <- function(tvec, dcsx, dcsy) {
  # s = int sqrt(dx/dt^2 + dy/dt^2) dt
  integral(function(t) sqrt(ppval(dcsx, t)^2 + ppval(dcsy, t)^2), tvec[1], tvec[2])
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
  # Plots
  # plot <- data.frame(x = seq(cs$breaks[1], cs$breaks[length(cs$breaks)], length = 100))
  # plot$y = ppval(cs, plot$x)
  # plot$cy = ppval(ccs, plot$x)
  # ggplot(plot, aes(x = x)) +
  #   geom_path(aes(y = y)) + geom_point(aes(y = y)) +
  #   geom_path(aes(y = cy)) + geom_point(aes(y = cy))
  return(ccs)
}

#--- Finite Difference Method ----
# offset = filter(dump$offset, onum == 1) %>% select(onum, nstep, norm, aveh, t_enum, offseth); print(offset)
FiniteDiff <- function(offset, var, order = NULL) {
  # https://en.wikipedia.org/wiki/Finite_difference_coefficient
  # Forward finite difference coefficients
  FFD1 <- list(
    c(     -1, 1                       ),
    c(   -3/2, 2, -1/2                 ),
    c(  -11/6, 3, -3/2,  1/3           ),
    c( -25/12, 4,   -3,  4/3, -1/4     ),
    c(-137/60, 5,   -5, 10/3, -5/4, 1/5)
  )
  # If no order specified, determine the order
  if (is.null(order)) {
    order = max(offset$nstep) - min(offset$nstep)
    if (order > 5) order = 5}
  # Group by offset number
  offset <- arrange(offset, onum, nstep)
  offset_list <- split(offset, offset$onum)
  # Determine gradient over each row
  offset_list <- lapply(offset_list, function(row) {
    if (sd(diff(row$norm)) > sqrt(.Machine$double.eps)) warning("Not equally sized steps")
    del = mean(diff(row$norm))
    diff = (sum(row[1:(order + 1), var] * FFD1[[order]]))/del
    return(cbind(
      row[1,], data.frame(diffvar = diff)))
  })
  # Clean up results
  offset_list <- bind_rows(offset_list) %>%
    select(onum, diffvar)
  colnames(offset_list)[ncol(offset_list)] <- paste0(var, "_diff")
  # Return the result
  return(LongJoin(offset, offset_list))
}

#--- Point in Polygon ----
PointinElement <- function(pts_df, poly_df, resize = NULL) {
  # Reduce the size of poly_df if possible
  if (!is.null(resize)) {
    poly_df <- filter(
      poly_df,
      x >= min(pts_df$x) - resize, y <= max(pts_df$x) + resize,
      y >= min(pts_df$y) - resize, y <= max(pts_df$y) + resize)}
  # Split the dataframe into a list based on enum and then remove enum from df in the list
  poly_list <- split(poly_df, poly_df$enum)
  # Convert the list to Polygon, then create a Polygons object
  poly_sp <- sapply(poly_list, function(poly){
    Polygons(list(Polygon(poly[, c("x", "y")])), ID = poly[1, "enum"])})
  poly_sp <- SpatialPolygons(poly_sp)
  # Convert points to coordinates
  pts_ps <- pts_df
  coordinates(pts_ps) <- ~x+y
  # Determine polygons points are in
  pts_return <- over(pts_ps, poly_sp, returnList = FALSE)
  pts_df$enum <- unique(poly_df$enum)[pts_return] 
  pts_df <- filter(pts_df, !is.na(enum))
  # Return
  return(pts_df)
}

