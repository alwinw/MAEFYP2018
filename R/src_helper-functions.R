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
  colnames(wallmesh) <- c("x", "y", "nxG", "nyG", "areaG", "elmt", "side")
  # out: wallmesh = data.frame(x, y, nxG, nyG, areaG)
  return(wallmesh)
}
# Determine the LE, TE, spline length, etc of the airfoil
# out:long_wall = data.frame(<wallmesh>, wnum, theta, up, wall, snum)
AirfoilLongWall <- function(wallmesh, 
                            dir = "clockwise") {
  # Note: Direction must be clockwise, cannot change!!
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
  if (dir != "clockwise")
    warning("Only clockwise supported for normals")
  # Patch LE and TE using elmt number
  # ASSUME that the LE corresponds to a node point!!
  long_wall <- long_wall %>%
    mutate(theta = 2*pi - theta,                                  # Reverse the direction for lower first
           up = theta >= pi) %>%
    group_by(elmt) %>%
    mutate(up = as.logical(round(mean(up))))                      # LE automatically patched
  long_wall$theta[long_wall$theta==2*pi & long_wall$up==FALSE] = 0
  long_wall <- long_wall %>% 
    mutate(avetheta = mean(theta)) %>%                            # Use average theta to "group" elements in sort by theta
    ungroup() %>%
    arrange(up, theta, avetheta) %>%
    select(-avetheta)
  # Output
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
    # cat("Determining spline distance\n")
    s <- cbind(data$s, lead(data$s))[-nrow(data),]                # Interval (s1, s2)
    s <- apply(s, 1, CubicSplineArcLength, dcsx, dcsy)
    # s <- pbapply(s, 1, CubicSplineArcLength, dcsx, dcsy)          # Arc length of interval
    s <- c(0, cumsum(s))                                          # Cumulative arc lenth
    error  = sum(data$s - s)                                      # Error
    data$s = s                                                    # Update loop
    # print(error)
  }
  csx <- cubicspline(data$s, data$x)                              # Cublic spline
  csy <- cubicspline(data$s, data$y)
  dcsx = CubicSplineCalc(csx, -1)                                 # Derivative of cubic splines
  dcsy = CubicSplineCalc(csy, -1)
  data$dxds <- ppval(dcsx, data$s)                                # Use cublic spline to determine derivatives
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
# Determine airfoil normals from spline
# out: long_wall = data.frame(<long_wall>, nxS, nyS, vecerr)
AirfoilNorm <- function(long_wall) {
  # Assume that nxG, nyG etc exist
  long_wall <- long_wall %>%
    mutate(nxS = nxG,
           nyS = -nxS/dydx) %>%
    mutate(distG = EucDist(nxG, nyG),
           distS = EucDist(nxS, nyS)) %>%
    mutate(nxG = nxG/distG,
           nyG = nyG/distG,
           nxS = nxS/distS,
           nyS = nyS/distS) %>%
    select(-distG, -distS) %>%
    mutate(vecerr = EucDist(nxG-nxS, nyG-nyS))
  # Checks
  if (min(cor(long_wall$nxG, long_wall$nxS), 
          cor(long_wall$nyS, long_wall$nyS)) < 0.999)
    warning("Low correlation between normals")
  if (sum(long_wall$vecerr > 0.005) > 0)
    warning(paste(sum(long_wall$vecerr > 0.005), "vector error > 0.5%"))
  # out: long_wall = data.frame(<long_wall>, nxS, nyS, vecerr)
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
    stringsAsFactors = FALSE,
    comment.char = "#")
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
#--- * Mass Data                                                  ----
# Load mass matrix
LoadMass <- function(seshpath,
                     extension = ".mass") {
  file           <- paste0(seshpath, extension)
  mass           <- read.table(file, skip = 0)
  colnames(mass) <- c("mass", "enum", "jnum")
  # out: mass = data.frame(x, y, enum, jnum)
  return(mass)
}
# Add mass to mesh
LongMass <- function(long_mesh, long_mass) {
  if (!all_equal(select(long_mesh,enum,jnum), 
                select(long_mass,enum,jnum)))
    warning("Unequal enum and jnum")
  long_mesh$mass <- long_mass$mass
  # out: mesh = data.frame(long_mesh, mass)
  return(long_mesh)
}
#--- * Wall Data                                                  ----
# Add enum and other important variables
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
#--- * Local Data                                                 ----
# Determine local mesh
# out: long_mesh = data.frame(<long_mesh>, local)
LocalMesh <- function(long_mesh, long_wall) {
  # First nodes and elements
  local_mesh <- filter(long_wall, node)
  local_mesh <- list(
    # nodes = unique(local_mesh$nnum),
    elements = unique(local_mesh$enum))
  local_mesh$mesh = long_mesh %>% 
    filter(node) %>% 
    select(enum, nnum) 
  local_mesh$mesh$local = 0
  local_mesh$mesh$local = 
    ifelse(local_mesh$mesh$enum %in% local_mesh$elements, 
           local_mesh$mesh$local - 1, 
           local_mesh$mesh$local)
  # Loop
  while (sum(local_mesh$mesh$local == 0) > 0) {
    local_mesh$nodes = unique(
      local_mesh$mesh[local_mesh$mesh$enum %in% local_mesh$elements,]$nnum)
    local_mesh$elements = unique(
      local_mesh$mesh[local_mesh$mesh$nnum %in% local_mesh$nodes,]$enum)
    local_mesh$mesh$local = 
      ifelse(local_mesh$mesh$enum %in% local_mesh$elements, 
             local_mesh$mesh$local - 1, 
             local_mesh$mesh$local)
  }
  local_mesh$mesh$local = local_mesh$mesh$local - min(local_mesh$mesh$local) + 1
  # Clean up and join
  local_mesh <- local_mesh$mesh %>% select(enum, local) %>% unique(.)
  long_mesh <- LongJoin(long_mesh, local_mesh)
  # out: long_mesh = data.frame(<long_mesh>, local)
  return(long_mesh)
}

#--- Dump File Calculation                                        ----
#--- * Dump Data                                                  ----
# Load GradFieldDump
# out: dump = list(time, kinvis, dump)
LoadGradFieldDump <- function(folder, dumpfile) {
  dumppath  <-  paste0(folder, dumpfile)
  filelines <- readLines(dumppath)
  # Read time
  time <-  filelines[grep("Time", filelines)[1]] %>%
    gsub("Time", "", .) %>%
    as.numeric(.)
  # Read Kinvis
  kinvis <-  filelines[grep("Kinvis", filelines)[1]] %>%
    gsub("Kinvis", "", .) %>%
    as.numeric(.)
  # Read table
  flowfield <- read.table(
    file = dumppath,
    skip = grep("ASCII", filelines),
    stringsAsFactors = FALSE)
  colnames(flowfield) <-                                          # o = omega, p = pressure
    c("u", "v", "p", "dodx", "dody", "dpdx", "dpdy", "o")
  # out: dump = list(time, kinvis, dump)
  return(list(time = time, kinvis = kinvis, dump = flowfield))
}
# Combine dump with mesh
# out: long_mesh = data.frame(<long_mesh>, <dump_dump>)
DumpMesh <- function(long_mesh, dump_dump) {
  # Check the number of rows are equal
  if (nrow(long_mesh) != nrow(dump_dump)) 
    warning("Unequal number of mesh numbers")
  # Join
  dump_dump$dnum <- 1:nrow(dump_dump)
  dump_dump <- cbind(long_mesh, dump_dump)
  # out: long_mesh = data.frame(<long_mesh>, <dump_dump>)
  return(dump_dump)
}
# Combine dump with wall
# out: dump_wall = data.frame(<long_wall>, <dump_dump>)
DumpWall <- function(long_wall, dump_dump) {
  dump_wall <- left_join(long_wall, select(dump_dump, -area),
                         by = c("x", "y", "enum", "node"))
  # out: dump_wall = data.frame(<long_wall>, <dump_dump>)
  return(dump_wall)
}

#--- * Accleration Data                                           ----
# Boundary equations from the session file
# out: none, function saved to environment
LoadSeshBCEq <- function(seshpath, bctext, bcfuncname = NULL) {
  # Read function from session file
  if (is.null(bcfuncname)) bcfuncname = bctext                    # Determine function name
  file = paste0(seshpath,".sesh")                                 # Session file name
  filelines <- readLines(file)
  bc <- filelines[grep(bctext, filelines)]
  # Clean up line
  bc %<>% gsub(bctext, "", .) %>%
    gsub("\t", "", .) %>%
    gsub("=", "", .) %>%
    gsub("PI", "pi", .)
  # Create the required bc function
  eval(parse(text = paste0(
    "BC_", tolower(bctext), " <- function(t) ", bc)),
    envir = .GlobalEnv)
}
# Relative Acceleration
# out: long_wall = data.frame(<long_wall>, as, an)
DumpAccel <- function(a, long_wall) {
  long_wall <- long_wall %>%
    mutate(
      a  = a,
      as = a*dxds/dydxlen,
      an = a*dyds/dydxlen)
  # out: long_wall = data.frame(<long_wall>, as, an)
  return(long_wall)
}
#--- * Pressure Data                                              ----
# out: dump_wall = data.frame(<dump_wall>, dpdsG, dpdsS)
DumpPres <- function(dump_wall, interp = TRUE) {
  # Wall Norm
  # Assumes clockwise direction
  dump_wall <- dump_wall %>%
    mutate(dpdsG = -nyG*dpdx +nxG*dpdy)
  if (interp) {
    # Cublic spline
    dump_pres <- dump_wall %>%
      select(x, y, s, p) %>%
      arrange(s) %>%
      unique(.)
    csp <- cubicspline(dump_pres$s, dump_pres$p)
    dcsp = CubicSplineCalc(csp, -1)
    dump_pres$dpdsS = ppval(dcsp, dump_pres$s)
    dump_wall <- LongJoin(dump_wall, dump_pres)
    # Checks
    if (cor(dump_wall$dpdsG, dump_wall$dpdsS) < 0.99)
      warning("Poor correlation between dpds methods")
    if (sum(abs(dump_wall$dpdsG - dump_wall$dpdsS)/dump_wall$dpdsG > 0.01))
      warning(paste(sum(abs(dump_wall$dpdsG - dump_wall$dpdsS)/dump_wall$dpdsG > 0.01)), " dpds error > 1%")
  }
  # out: dump_wall = data.frame(<dump_wall>, dpdsG, dpdsS)
  return(dump_wall)
}
#--- * Vorticity Data                                             ----
# Determine normal offsets
# out: dump_offs = tibble(x, y, nxS, nyS, aveh, wnum, onum, nstep, offseth, norm, wall)
AirfoilOffset <- function(dump_wall, 
                          totdist = 0.008, nsteps = 5, varh = TRUE, scale = 1) {
  # Variables of interest
  dump_offs <- dump_wall %>%
    select(x, y, theta, s, nxG, nyG, nxS, nyS, dodx, dody, aveh, enum, wnum) %>%
    group_by(x, y, theta) %>%
    mutate(aveh = mean(aveh)) %>%
    ungroup() %>%
    select(-theta)
  dump_offs$onum <- 1:nrow(dump_offs)
  # Repeat for nsteps
  dump_offs <- slice(dump_offs, rep(1:n(), each = nsteps + 1))
  dump_offs$nstep = rep(0:nsteps, length.out = nrow(dump_offs))
  dump_offs <- dump_offs %>%
    mutate(
      offseth = ifelse(is.na(aveh) | !varh, totdist, aveh)*scale,
      norm = offseth*nstep/nsteps,
      x = x - nxS*norm,
      y = y - nyS*norm) %>%
    mutate(wall = ifelse(nstep == 0, TRUE, FALSE),
           wnum = ifelse(wall, wnum, NA))
  # out: dump_offs = tibble(x, y, nxS, nyS, aveh, wnum, onum, nstep, offseth, norm, wall)
  return(dump_offs)
}
# Interpolate vorticity data onto offset points
# dump_offs = dump$off; dump_dump = dump$dump; localmax = 2;

DumpVortInterp <- function(dump_offs, dump_dump,
                           localmax = 2, linear = TRUE, extrap = FALSE, round = NULL) {
  dump_offs$enumo <- dump_offs$enum
  # Points (ps) in Spacial Polygon (ps)
  poly_df <- dump_dump %>%
    filter(local <= localmax, node) %>%
    arrange(enum, ncorner) %>%
    select(x, y, enum)
  poly_sp <- split(poly_df, poly_df$enum)
  poly_sp <- sapply(poly_sp, function(poly) {
    Polygons(list(Polygon(poly[, c("x", "y")])), ID = poly[1, "enum"])})
  poly_sp <- SpatialPolygons(poly_sp)
  pts_ps  <- dump_offs %>%
    select(x, y)
  coordinates(pts_ps) <- ~x + y
  pts_ret <- over(pts_ps, poly_sp, returnList = FALSE)
  dump_offs$enum <- unique(poly_df$enum)[pts_ret] 
  if (sum(is.na(dump_offs$enum)) > 0) {
    warning("Some points not found in polygon, maybe larger localmax?")
    warning("Patching based on original mesh number...")
    dump_offs$enum <- ifelse(is.na(dump_offs$enum), dump_offs$enumo, dump_offs$enum)
  }
  # Interpolation
  mesh <- dump_dump %>%
    filter(local <= localmax) %>%
    select(x, y, o, enum)
  mesh_list <- split(mesh, mesh$enum)
  offs <- dump_offs %>%
    select(x, y, enum) %>%
    arrange(enum) %>%
    as.data.frame(.)
  offs_list <- split(offs, offs$enum)
  # SPECIAL NOTE: since enum is int, mesh_list[[57]] =/= mesh_list[["57"]]
  # mesh_list[[names(offs_list)[1]]] =/= mesh_list[[57]]
  # Testing
  if (FALSE) {
    mesh_val <- mesh_list[[names(offs_list)[1]]]
    mesh_val <- round(mesh_val, 8)
    offs_val <- offs_list[[names(offs_list)[1]]]
    offs_val <- mesh_val[7,] #offs_val[12,]
    ggplot(mapping = aes(x, y)) + 
      geom_point(shape = 19, colour = "red", size = 3, alpha = 0.5,
                 data = mesh_val[7,]) + 
      geom_point(aes(colour = o), shape = 19, data = mesh_val) +
      geom_point(shape = 21, data = offs_val) +
      geom_polygon(fill = NA, colour = "black", 
                   data = mesh_val[chull(mesh_val[,1], mesh_val[,2]),])
    # Interpolate(mesh_val, offs_val, TRUE)
    interpp( x = mesh_val[,1],  y = mesh_val[,2],  z = mesh_val[,3],
            xo = offs_val[,1], yo = offs_val[,2],
            linear = FALSE,
            extrap = FALSE,
            duplicate = "strip")
    }
  intp_out  <- lapply(names(offs_list), function(enum) {
    intp <- Interpolate(mesh_list[[enum]], offs_list[[enum]], 
                        linear = linear, extrap = extrap, round = round)
  })
  intp_out <- bind_rows(intp_out)
  if (sum(abs(intp_out$x-offs$x) + abs(intp_out$y-offs$y)))
    warning("enum is out of order")
  offs <- cbind(intp_out, enum=offs$enum)
  # Combine back
  dump_offs <- LongJoin(dump_offs, offs)
  # out: dump_offs = tibble(<dump_offs>, enumo, o)
  return(dump_offs)
}
# Create the vorticity gradients for both methods
DumpVortGrad <- function(dump_offs) {
  dump_offs <- dump_offs %>%
    rename(dodzS = o_diff) %>%
    mutate(dodzG = -nxG*dodx - nyG*dody)
  # BROKEN, THERE ARE NA FROM THE INTERPOLATION!!
  check <- dump_offs %>% 
    filter(!is.na(dodzS)) %>%
    select(dodzG, dodzS) %>%
    unique(.)
  # print(cor(check$dodzG, check$dodzS))
  # print(max(abs(check$dodzG - check$dodzS)/check$dodzG))
  # print(median(abs(check$dodzG - check$dodzS)/check$dodzG))
  return(dump_offs)
}

DumpVortJoin <- function(dump_wall, dump_offs, dump_kinvis) {
  # Combine offset data into dump_wall
  dump_offs <- filter(dump_offs, nstep==0)
  if (sum(abs(dump_wall$s - dump_offs$s)) > 0)
    warning("Difference in s for wall and offs detected")
  dump_offs <- rename(dump_offs, interpo = o)
  dump_wall <- cbind(
    dump_wall,
    select(dump_offs, offseth, enumo, interpo, dodzS, dodzG))
  dump_wall <- dump_wall %>% 
    mutate(LHSS   = -as + dpdsS,
           RHSS   = -dodzS*dump_kinvis,
           nserrS =  LHSS - RHSS,
           LHSG   = -as + dpdsG,
           RHSG   = -dodzG*dump_kinvis,
           nserrG =  LHSG - RHSG)
  return(dump_wall)
}

DumpVortOnly <- function(dump_wall, dump_kinvis) {
  dump_wall <- dump_wall %>% 
    mutate(dodzG = -nxG*dodx - nyG*dody) %>% 
    mutate(LHSG   = -as + dpdsG,
           RHSG   = -dodzG*dump_kinvis,
           nserrG =  LHSG - RHSG)
}

DumpFlow <- function(dump_dump) {
  dump_infl <- dump_dump %>% 
    filter(x == min(dump_dump$x)) %>% 
    select(u, v)
  dump_outf <- dump_dump %>% 
    filter(x == max(dump_dump$x)) %>% 
    select(u, v)
  
}

#--- * Integral                                                   ----
LoadIntegral <- function(folder, dumpfile, vars = c("u", "v","o")) {
  # Ready the file into filelines
  dumppath  <-  paste0(folder, strsplit(dumpfile, ".dump")[[1]], ".integral")
  filelines <- readLines(dumppath)
  # Manipulate to get a useful data.frame
  varlines <- c(1, grep("centroid", filelines))
  vartable <- data.frame(
    start = varlines[-length(varlines)],
    end   = varlines[-1] )
  varlines <- data.frame(char = filelines[varlines[-1]]) 
  vardf <- varlines%>% 
    separate(char, c("var", "b", "c"), ":") %>% 
    separate(b, c("intvar", "junk"), ",")   %>% 
    separate(c, c("centx", "centy"), ",")
  vardf$intvar <- as.numeric(vardf$intvar)
  vardf$centx  <- as.numeric(vardf$centx)
  vardf$centy  <- as.numeric(vardf$centy)
  vardf <- cbind(vardf, vartable)
  rm(varlines, vartable)
  # Read each variable of interest (vars vector)
  vardf <- filter(vardf, var %in% vars)
  varli <- split(vardf, vardf$var)
  varoutp <- lapply(varli, function(varline) {
    values = read.table(file = dumppath, 
               skip = varline$start, nrows = varline$end - varline$start - 1)
    values = unique(values)
    # colnames(values) <- c("enum", paste0("int.", varline$var))
    values = values[,2]
    return(values)
  })
  varoutp <- bind_cols(varoutp) %>% 
    mutate(enum = row_number())
  vardf <- select(vardf, -start, -end, -junk)
  
  # Return output
  return(list(total = vardf, elem = varoutp))
}

DumpIntegral <- function(dump_dump, long_enum) {
  # Variables of interest
  inte_dump <- dump$dump %>% 
    select(x, y, o) %>% 
    mutate(ox = o*(x-3), oy = -o*y) %>% 
    select(-x, -y)
  # Multiply by mass
  inte_dump <- inte_dump*dump_dump$mass
  inte_dump$enum = dump_dump$enum
  # Integrate per element (cannot do per point since not really "dA")
  inte_dump <- inte_dump %>% 
    group_by(enum) %>% 
    summarise_all(.funs = "sum")
  # Total integration
  inte_totl <- inte_dump %>% 
    select(-enum) %>% 
    summarise_all(.funs = "sum")
  
}

#--- Numerical Methods ----
# Heavyside function (step function)
heav <- function(t) ifelse(t>0,1,0)

# Distance function
LagDist <- function(x, y) sqrt((x - lag(x))^2 + (y - lag(y))^2)
EucDist <- function(x, y) sqrt(x^2 + y^2)

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
CubicSplineArcLength <- function(tvec, dcsx, dcsy, baseint = TRUE) {
  # s = int sqrt(dx/dt^2 + dy/dt^2) dt
  if (baseint) {
    # Use base integration method (much faster based in benchmarking!)
    integrate(function(t) sqrt(ppval(dcsx, t)^2 + ppval(dcsy, t)^2), tvec[1], tvec[2])$value
  } else {
    # Use pracma integration (must slower)
    integral(function(t) sqrt(ppval(dcsx, t)^2 + ppval(dcsy, t)^2), tvec[1], tvec[2])
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

#--- Interpolation ----
Interpolate <- function(mesh, input, 
                        linear = FALSE, extrap = FALSE, round = NULL) {
  names <- colnames(mesh)
  if (!is.null(round))
    mesh <- round(mesh, round)
  # Interpolate
  suppressWarnings(
    output <- as.data.frame(
      interpp( mesh[,1],  mesh[,2],  mesh[,3],
              input[,1], input[,2],
              linear = linear,
              extrap = extrap,
              duplicate = "strip"))
  )
  colnames(output) <- names[1:3]
  # Check for extremely large/small values and NA
  output[,3] = ifelse(output[,3] < min(mesh[,3]) | output[,3] > max(mesh[,3]),
                      NA, output[,3])
  if (sum(is.na(output[,3])) > 0) warning("NA found in bicubic interpolation")
  return(output)
}

#--- Batch                                                        ----
ListSesh <- function(batchfolder, pathsplit = c("airfoil", "seshname")) {
  # List of files in folder (should be made more robust later)
  batchlist <- list.files(batchfolder, pattern = ".sesh", recursive = TRUE)
  # Convert character to dataframe
  batchlist <- data.frame(path = unlist(strsplit(batchlist, "*.sesh"))) %>%
    separate(path, pathsplit, sep = "/")
  if (length(pathsplit) != 1) {
    # Split the folder
    folder <- batchlist[colnames(batchlist)[colnames(batchlist) != "seshname"]] %>% 
      unite(folder, sep = "/")
    batchlist$folder <- folder$folder
    # Create output
    batchlist <- batchlist %>% 
      mutate(folder = paste0(batchfolder, "/", folder, "/")) %>%
      mutate(seshpath = paste0(folder, seshname)) %>% 
      select(airfoil, seshname, folder, seshpath)
  } else {
    # Create output
    batchlist <- batchlist %>% 
      mutate(folder = paste0(batchfolder, "/")) %>%
      mutate(seshname = airfoil) %>% 
      mutate(seshpath = paste0(folder, seshname)) %>% 
      select(airfoil, seshname, folder, seshpath)
  }
  # Return list of session files
  return(batchlist)
}

# List dump files in directory
ListDump <- function(folder, seshname) {
  # List of files in folder (should be made more robust later)
  dumplist <- list.files(folder, pattern = paste0(seshname,"-"))
  # Return file list
  return(dumplist)  # character
}

LoadSeshTokenWords <- function(seshpath, tokenwords) {
  # Session file name
  file = paste0(seshpath,".sesh")
  # Read the session file to grep lines later
  filelines <- readLines(file)
  # For each keyword load the associated data
  output <- lapply(tokenwords, function(tokenword) {
    value = filelines[grep(tokenword, filelines)] %>%
      gsub(tokenword, "", .) %>%
      gsub("\t", "", .) %>%
      gsub("=", "", .) %>%
      as.numeric(.)
    return(data.frame(tokenword = tokenword,
                      tokenvalue = value,
                      stringsAsFactors = FALSE))
  })
  # Convert to data.frame
  output <- bind_rows(output)
  # Return the session file contents
  return(output)   # List
}

#--- Plot Outputs                                                 ----
#--- * Plot Setup                                                 ----
# Variables used in multiple plots
PlotSetup <- function(plot_wall, plot_data) {
  # Vertical lines to make TE, LE, TE
  plot_vlines <- data.frame(
    telo = min(plot_wall$s),
    le   = plot_wall[plot_wall$theta == pi,]$s[1],
    teup = max(plot_wall$s))
  # Location of surface labels
  plot_surf <- data.frame(
    x = c(plot_vlines$telo,
          mean(c(plot_vlines$telo, plot_vlines$le)),
          plot_vlines$le,
          mean(c(plot_vlines$le, plot_vlines$teup)),
          plot_vlines$teup),
    y = rep(-40, 5),
    labels = c("TE", "Lower", "LE", "Upper", "TE"))
  # xbreak locations on x-axis
  plot_xbreaks <-
    c(seq(plot_vlines$telo, plot_vlines$le, length.out = 3)[1:2],
      seq(plot_vlines$le, plot_vlines$teup, length.out = 3))
  # Title of plot with key information
  plot_title <- paste0(
    plot_data$airfoil, "\n",
    paste("Time:",                sprintf("%05.3f",  plot_data$time  )), "   ",
    paste("Acceleration:",        sprintf("%+07.4f", plot_data$a   )), "\n",
    paste("Kinematic Viscosity:", sprintf("%.4f",   plot_data$kinvis))
  )
  # Plot filename
  plot_filename <- paste(
    plot_data$ID,
    paste0("v", format(plot_data$kinvis, scientific = TRUE)),
    paste0("t", sprintf("%05.3f",  plot_data$time)),
    paste0("a", sprintf("%+07.3f", plot_data$a)),
    # "ns",
    sep = "_")
  # Plots limits
  # Add limits for airfoil, LE, TE, etc here
  # Output
  plot_setup <- list(
    vlines   = plot_vlines,
    surf     = plot_surf,
    xbreaks  = plot_xbreaks,
    title    = plot_title,
    filename = plot_filename)
  return(plot_setup)
}
