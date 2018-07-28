#============================#
# Analysis of NACA0012 Example
# Alwin Wang
#----------------------------#

#--- Set Up                                                       ----
# Use rstudioapi to get saved location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
#--- ^ Scripts                                                    ----
# Source Required Scripts
source("src_library-manager.R")                                 # Call libraries and install missing ones
source("src_numerical-methods.R")                               # Load custom numerical methods
# Additional scripts here
# ggplot2 setup (consider moving to a separate script)
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 
#--- ^ * Required Paths                                           ----
saveplot   = "../src-example/NACA0012/"
airfoil    = "NACA0012-AoA04"
folderpath = "../src-example/NACA0012/results/"
seshpath   = "RE-10000-sine-0.001-2000"
dumppath   = "RE-10000-sine-0.001-2000-03.dump"
# Required input dataframes
airfoildata <- data.frame(
  airfoil  = airfoil, 
  seshname = seshpath,
  folder   = folderpath,
  seshpath = paste0(folderpath, seshpath),
  stringsAsFactors = FALSE)
meshdata <- data.frame(
  airfoildata,
  tokenword  = "N_P",
  tokenvalue = 5,
  ID         = "NACA0012-AoA04-N_P5",
  stringsAsFactors = FALSE)
dumpdata <- data.frame(
  meshdata,
  dumpfile = dumppath,
  dumppath = paste0(folderpath, dumppath),
  stringsAsFactors = FALSE)
rm(airfoil, folderpath, seshpath, dumppath)

#--- Airfoil Calculation                                          ----
# Determine things like spline distance once per unique airfoil
# i.e. only load the boundary profile and wall output
#--- ^ Boundary Data                                              ----
#--- ^ * LoadBndry(airfoildata)                                   ----
# out:  bndry = data.frame(x, y)
file = paste0(airfoildata$folder, "bndry_prf", ".dat")
bndry           <- read.table(file)
colnames(bndry) <- c("x", "y")
  rm(file)
  # ggplot(bndry, aes(x, y)) + geom_point() + coord_fixed()
#--- ^ Wall Mesh Data                                             ----
#--- ^ * LoadWallGrad(airfoildata)                                ----
# out:  wallmesh = data.frame(x, y, nxG, nyG, areaG)
file = paste0(airfoildata$seshpath, ".wallgrad")
wallmesh           <- read.table(file, skip = 1)
colnames(wallmesh) <- c("x", "y", "nxG", "nyG", "areaG")
  rm(file)
  # wallmeshplot <- rbind(wallmesh[,1:2], wallmesh[,1:2] - wallmesh[,3:4])
  # wallmeshplot$wnum <- 1:nrow(wallmesh)
  # ggplot(wallmeshplot, aes(x, y, group = wnum)) + geom_line()
  #   rm(wallmeshplot)
#--- ^ * AirfoilLongWall(wallmesh)                                ----
# out:  long_wall = data.frame(<wallmesh>, wnum, theta, up, wall, snum)
te = wallmesh[which.max(wallmesh$x), 1:2]                       # Most far right point
le <- wallmesh %>%                                              # Furthest point from te
  mutate(dist = sqrt((x-te$x)^2 + (y-te$y)^2)) %>% arrange(-dist)
le = le[1, 1:2]
cp = (le+te)/2
tetheta = atan2(te$y-cp$y, te$x-cp$x)
wallmesh$wnum = 1:nrow(wallmesh)
long_wall <- wallmesh %>%
  mutate(theta = atan2(y - cp$y, x - cp$x) - tetheta) %>%
  mutate(theta = theta + ifelse(theta < 0, 2*pi, 0)) %>%        # CANNOT use sign(theta) since theta can be zero
  arrange(theta, wnum)
# Patch TE 
long_wall$theta[1] = long_wall$theta[1] + 2*pi
long_wall <- long_wall  %>%
  mutate(theta = 2*pi - theta) %>%                              # TE -> lower -> LE -> upper -> TE
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
long_wall <- mutate(long_wall, wall = !is.na(wnum))             # Add helpful variables
long_wall$snum = 1:nrow(long_wall)
  rm(cp, le, lepatch, te, tetheta)
  # ggplot(long_wall, aes(x, y, colour = theta, shape = up)) + geom_point() + geom_path() + coord_fixed()
#--- ^ * AirfoilSpline(long_wall)                                 ----
# out: long_wall = data.frame(<long_wall>, s, dxds, dyds, dydx, dydxlen)
data         = long_wall[,colnames(long_wall) %in% c("x", "y", "theta")]
originalnrow = nrow(long_wall)
# Note: spline is closed by repeated the first point, need a robust method for TE
data = unique(data)                                             # Remove duplicate rows
# Iterate to find the best spline distance
data$s = c(0, LagDist(data$x,data$y)[-1])
data$s = cumsum(data$s)
error  = 1                                                      # Initialise variable
while (abs(error) > 0.001) {
  csx <- cubicspline(data$s, data$x)                            # Cublic spline
  csy <- cubicspline(data$s, data$y)
  dcsx = CubicSplineCalc(csx, -1)                               # Derivative of cubic splines
  dcsy = CubicSplineCalc(csy, -1)
  cat("Determining spline distance\n")
  s <- cbind(data$s, lead(data$s))[-nrow(data),]                # Interval (s1, s2)
  s <- pbapply(s, 1, arclength)                                 # Arc length of interval
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
data$dyds    <- ppval(dcsy, data$s)
data$dydx    <- data$dyds/data$dxds
data$dydxlen <- sqrt(data$dxds^2 + data$dyds^2)
# Combine data back with long_wall
long_wall <- full_join(long_wall, data, by = c("x", "y", "theta"))
if (originalnrow != nrow(long_wall)) {
  warning("Merge Failed")}                                      # Check the number of rows is correct
rm(data, csx, csy, dcsx, dcsy, error, originalnrow, s)
#--- > Airfoil Calc Output                                        ----
airfoildatalist <- list(
  airfoil   = airfoildata,
  bndry     = bndry,
  long_wall = long_wall)
rm(airfoildata, bndry, long_wall)
#--- Session and Mesh Calculation                                 ----
#--- ^ Airfoil Data                                               ----
long <- list()
long$walldata = airfoildatalist$long_wall
#--- ^ Session Data                                               ----
keywords <- list(                                               # Keywords in session to read
  c("NODES", "nnum", "x", "y", "z"),
  c("ELEMENTS", "enum", "shapetag", "n1", "n2", "n3", "n4", "junk"),
  c("SURFACES", "snum", "element", "side", "bctag", "bc", "junk"),
  c("CURVES", "cnum", "element", "side", "curvetag", "curvedata", "junk"))
#--- ^ * LoadSeshFileKeywords(meshdata, keywords)                 ----

