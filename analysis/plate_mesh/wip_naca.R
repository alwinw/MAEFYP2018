#============================#
# Mesh Generation
# Alwin Wang
#----------------------------#

# LOWER SURFACE IS NOT CORRECT FOR CAMBERED AIRFOILS

#--- Set Up ----
# Use rstudioapi to get saved location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
#--- * Scripts ----
# Source Required Scripts
srcpath = "../../R/"
source(paste0(srcpath, "src_library-manager.R"))                # Call libraries and install missing ones
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 

#--- Inputs ----
input <- list(
  x = data.frame(min = -1.75, max = 5.4375, n = 61, polyn = 4), #51
  y = data.frame(min = -1.25, max = 1.25,   n = 21, polyn = 3), #21
  t = 0.12,
  alpha = 4, # degrees
  NACA = 0012
)

#--- Polynomial ----
c = data.frame(x = 0.5*cos(input$alpha*pi/180),
               y = 0.5*sin(input$alpha*pi/180))
#--- * x ----
dx = (input$x$max - input$x$min)/input$x$n
listx = list(
  LE = seq(input$x$min, -c$x, length.out = round((-c$x - input$x$min)/dx) + 1),
  CD = seq(-c$x, c$x, length.out = 2*round(c$x*2/dx) + 1),
  TE = seq(c$x, -input$x$min, length.out = round((-input$x$min - c$x)/dx) + 1) )
polyx = list(
  LE = polypow(c(1,  c$x), input$x$polyn + 1), 
  CD = polyint(polymul(polypow(c(1,-c$x),input$x$polyn), polypow(c(1,c$x),input$x$polyn))),
  TE = polypow(c(1, -c$x), input$x$polyn + 1) )
polyx$LE = polyadd(polyx$LE/((input$x$min + c$x)^input$x$polyn), c(-c$x))
polyx$CD = c$x/(polyval(polyx$CD, c$x)) * polyx$CD
polyx$TE = polyadd(polyx$TE/((-input$x$min - c$x)^input$x$polyn), c( c$x))
outpx = list(
  LE = polyval(polyx$LE, listx$LE),
  CD = polyval(polyx$CD, listx$CD),
  TE = polyval(polyx$TE, listx$TE) )
outpx$WK = seq(outpx$TE[length(outpx$TE)], input$x$max, 
               length.out = round((input$x$max - outpx$TE[length(outpx$TE)])/
                                    (outpx$TE[length(outpx$TE)] - outpx$TE[length(outpx$TE)-1])))
#--- * y ----
listy = list(
  LO = seq(input$y$min, c$y, length.out = round(input$y$n/2) + 1),
  UP = seq(c$y, input$y$max, length.out = round(input$y$n/2) + 1) )

polyy = list(
  LO = polypow(c(1, -c$y), input$y$polyn + 1), 
  UP = polypow(c(1, -c$y), input$y$polyn + 1) )
polyy$LO = polyadd(polyy$LO/((input$y$max + c$y)^input$y$polyn)*(-2*mod(input$y$polyn, 2) + 1), c(c$y))
polyy$UP = polyadd(polyy$UP/((input$y$max - c$y)^input$y$polyn), c(c$y))
outpy = list(
  LO = polyval(polyy$LO, listy$LO),
  UP = polyval(polyy$UP, listy$UP) )
#--- * Plot ----
if (FALSE) {
  ggplot() +
    geom_vline(xintercept = outpx$LE, colour = "blue" ) +
    geom_vline(xintercept = outpx$TE, colour = "blue" ) +
    geom_vline(xintercept = outpx$CD, colour = "red"  ) +
    geom_vline(xintercept = outpx$WK, colour = "green") +
    geom_hline(yintercept = outpy$LO, colour = "red"  ) +
    geom_hline(yintercept = outpy$UP, colour = "blue" ) +
    coord_fixed(xlim = c(input$x$min, input$x$max),
                ylim = c(input$y$min, input$y$max), 
                expand = FALSE)
}


#--- Output ----
#--- * Lists ----
# x
dfx <- data.frame(
  x   = c(outpx$LE, outpx$CD[-1], outpx$TE[-1], outpx$WK[-1]),
  bcx = c("v" , rep("FD", length(outpx$LE) - 2), 
          "LE", rep("CD", length(outpx$CD) - 2),
          "TE", rep("BD", length(outpx$TE) - 2),
          rep("BD", length(outpx$WK) - 1),
          "o") )
# y
dfy <- data.frame(
  y   = c(outpy$LO, outpy$UP[-1]),
  bcy = c("w" , rep("LO", length(outpy$LO) - 2),
          "CD", rep("UP", length(outpy$UP) - 2),
          "s") )
# y' for trailing edge
dfy = dfy[nrow(dfy):1,]
dfyd <- data.frame(
  y   = -rev(dfy$y),
  bcy = dfy$bcy )

# Airfoil
AirfoilCurve <- function(NACA, x = 0, a = -0.5, c = 1, out = "all") {
  # Max camber; Location of m; Thickness
  m = (NACA %/% 1000) / 100
  p = (NACA %/% 100 %% 10) / 10
  t = (NACA %% 100) / 100
  # Test if x is within range
  on = ifelse(x >= a & x <= (a + c), 1, 0) # allows for root-finding
  # Determine the camber line yc
  yc = ifelse(x < p * c + a, 
              m/p^2 * (2*p*((x-a)/c) - ((x-a)/c)^2),
              m/(1-p)^2  * (1 - 2*p + 2*p*((x-a)/c) - ((x-a)/c)^2)
  )
  # Determine the gradient of the camber line dycdx
  dycdx = ifelse(x < p * c + a,
                 2*m/p^2 * (p - (x-a)/c),
                 2*m/(1-p)^2 * (p - (x-a)/c)
  )
  # Determine the magnitude and direction of the thickness
  theta = atan(dycdx)
  yt = 5*t*(0.2969*sqrt(abs((x-a)/c)) - 0.1260*((x-a)/c) - 0.3516*((x-a)/c)^2 +
              0.2843*((x-a)/c)^3 - 0.1036*((x-a)/c)^4)
  # Add the thickness to the camber line
  xU = x - yt*sin(theta)
  yU = yc + yt*cos(theta)
  xL = x + yt*sin(theta)
  yL = yc - yt*cos(theta)
  # Output depending on the Out parameter
  if(out == "all")
    summary = data.frame(x, yc, dycdx, theta, yt,  xU, yU,  xL, yL)
  else if(out == "coord")
    summary = data.frame(x, xU, yU, xL, yL)
  else if (out == "upper")
    summary = data.frame(x = xU, y = yU)
  else if (out == "lower")
    summary = data.frame(x = xL, y = yL)
  # Return the output
  return(summary * on)
}

Airfoilx <- function(NACA, xO, surf = "upper", a = -0.5, c = 1, tol = 1e-9, out = "x") {
  # Use the rooting finding in {stats} to find the root
  rootfind <- uniroot(function(x) AirfoilCurve(NACA, x, a, c, out = surf)$x - xO,
                      lower = a, upper = a + c,
                      tol = tol)
  if(out == "x")
    return(rootfind$root)
  else if(out ==  "all")
    return(rootfind)
  else if(out == "str")
    return(str(rootfind))
}

Airfoilxd <- function(NACA, alpha, x0, surf = "upper", a = -0.5, c = 1, tol = 1e-9, out = "x") {
  rootfind <- uniroot(
    function(x) {
      local  = AirfoilCurve(NACA, x, a, c, out = surf)
      r.x    = local$x*cos(alpha) + local$y*sin(alpha)
      return(r.x - x0) },
    lower = a, upper = a + c,
    tol = tol
  )
  if(out == "x")
    return(rootfind$root)
  else if(out ==  "all")
    return(rootfind)
  else if(out == "str")
    return(str(rootfind))
}

xa = -0.5; xc = 1
xvec = seq(0.5, -0.5, -0.025)
xvec = -2*xa/xc^3 * (xvec - xa)^3 + xa
airfoil <- data.frame(x = xvec)
airfoil$x0U <- sapply(airfoil$x, function(x) Airfoilx(input$NACA, x, "upper"))
airfoil$x0L <- sapply(airfoil$x, function(x) Airfoilx(input$NACA, x, "lower"))
airfoil$yU  <- AirfoilCurve(input$NACA, airfoil$x0U)$yU
airfoil$yL  <- AirfoilCurve(input$NACA, airfoil$x0U)$yL
airfoil <- data.frame(
  x0 = c(airfoil$x, rev(airfoil$x)[-1]),
  y0 = c(airfoil$yU, rev(airfoil$yL)[-1])
)
airfoil <- airfoil %>%  
  mutate(
    x =   x0*cos(input$alpha*pi/180) + y0*sin(input$alpha*pi/180),
    y = - x0*sin(input$alpha*pi/180) + y0*cos(input$alpha*pi/180) ) %>% 
  select(-x0, -y0)

if(FALSE) {
  ggplot(airfoil, aes(x, y)) +
    geom_path() +
    coord_fixed()
}

# y for x where x in CD
flags <- function(n1, i, j) {
  if (n1$bcx == "CD") {
    # Linearly interpolate to change
    n1$y = (dfyd$y[j] - n1$y)/(c$x*2)*(n1$x + c$x) + n1$y
  } else if (n1$bcx %in% c("TE", "BD", "o")) {
    n1$y = dfyd$y[j]
  }
  return(n1)
}

# airfoil
AirfoilTrans <- function(n1) {
  if (n1$bcx == "CD") {
    # Add on airfoil thickness
    alpha = input$alpha*pi/180
    x     = n1$x             # x location of current node
    y     = n1$y             # y location of current node
    y0    = -n1$x*tan(alpha) # y location of chord line
    ym    = input$y$max*0.8  # maximum y value
    if (n1$bcy == "UP" | n1$bcy == "CD") {
      x0 = Airfoilxd(input$NACA, alpha, x, "upper")      # x value on chord line for x value on surface
      ys = AirfoilCurve(input$NACA, x0)$yU/cos(alpha)
      dy = ys*(cos(pi/(2*ym)*(y-y0)))^2
      n1$y = n1$y + dy
    }
    if (n1$bcy == "LO") {
      x0 = Airfoilxd(input$NACA, alpha, x, "lower")
      ys = AirfoilCurve(input$NACA, x0)$yL/cos(alpha)
      dy = ys*(cos(pi/(2*ym)*(y-y0)))^2
      n1$y = n1$y + dy
    }
  }
  return(n1)
}

#--- * Nodes ----
# Numbering
nx = nrow(dfx)
ny = nrow(dfy)
matnnum <- matrix(1:(nx*ny), nrow = ny, ncol = nx, byrow = TRUE)
dfnode <- data.frame()
for (j in 1:ny) {
  for (i in 1:nx) {
    dfrow = data.frame(
      nnum = matnnum[j, i], x = dfx$x[i], y = dfy$y[j], z = 0, bcx = dfx$bcx[i], bcy = dfy$bcy[j])
    dfrow = flags(dfrow, i, j)
    dfrow = AirfoilTrans(dfrow)
    dfnode = rbind(dfnode, dfrow)
  } 
}
#--- * Plot ----
if (FALSE) {
  ggplot(dfnode, aes(x, y, colour = nnum)) +
    geom_point() + 
    scale_colour_gradientn(colours = spectralpalette(10))
}

#--- * Elements -----
enum = 0
dfelem <- data.frame()
surf = 0
dfsurf <- data.frame()

for (j in 1:(ny - 1)) {
  for (i in 1:(nx - 1)) {
    # 1 ---- 4
    # | enum |
    # 2 ---- 3
    enum = enum + 1
    # Determine node corners and boundary conditions
    n1 = data.frame(x = dfx$x[i  ], y = dfy$y[j  ], nnum = matnnum[j  , i  ], bcx = dfx$bcx[i  ], bcy = dfy$bcy[j  ])
    n2 = data.frame(x = dfx$x[i  ], y = dfy$y[j+1], nnum = matnnum[j+1, i  ], bcx = dfx$bcx[i  ], bcy = dfy$bcy[j+1])
    n3 = data.frame(x = dfx$x[i+1], y = dfy$y[j+1], nnum = matnnum[j+1, i+1], bcx = dfx$bcx[i+1], bcy = dfy$bcy[j+1])
    n4 = data.frame(x = dfx$x[i+1], y = dfy$y[j  ], nnum = matnnum[j  , i+1], bcx = dfx$bcx[i+1], bcy = dfy$bcy[j  ])
    # Determine wall bcs
    bc1 = NA; bc2 = NA; bc3 = NA; bc4 = NA;
    if (j == 1) {         # TOP | side 4
      bc4 = "s"; surf = surf + 1
      dfsurf = rbind(dfsurf, data.frame(
        i = i, j = j, surf = surf, enum = enum, side = 4, bc = bc4, nnum1 = n4$nnum, nnum2 = n1$nnum)) }
    if (i == (nx - 1)) {  # RIGHT | side 3
      bc3 = "o"; surf = surf + 1
      dfsurf = rbind(dfsurf, data.frame(
        i = i, j = j, surf = surf, enum = enum, side = 3, bc = bc3, nnum1 = n3$nnum, nnum2 = n4$nnum)) }
    if (j == (ny - 1)) {  # BOTTOM | side 2
      bc2 = "w"; surf = surf + 1
      dfsurf = rbind(dfsurf, data.frame(
        i = i, j = j, surf = surf, enum = enum, side = 2, bc = bc2, nnum1 = n2$nnum, nnum2 = n3$nnum)) }
    if (i == 1) {         # LEFT | side 1
      bc1 = "v"; surf = surf + 1
      dfsurf = rbind(dfsurf, data.frame(
        i = i, j = j, surf = surf, enum = enum, side = 1, bc = bc1, nnum1 = n1$nnum, nnum2 = n2$nnum)) }
    
    # Side 2 airfoil bc
    if (n2$bcy == "CD" & n3$bcy == "CD") { # Element above plate | side 2
      if (n3$bcx %in% c("CD", "TE")) {
        # if (n2$bcx == "CD") {n2$bcy = "UP"}
        bc2 = "p"; surf = surf + 1
        dfsurf = rbind(dfsurf, data.frame(
          i = i, j = j, surf = surf, enum = enum, side = 2, bc = bc2, nnum1 = n2$nnum, nnum2 = n3$nnum)) }
    }
    # Side 4 airfoil bc
    if (n4$bcy == "CD" & n1$bcy == "CD") {
      if (n4$bcx %in% c("CD", "TE")) { # Element below plote | side 4
        # if (n1$bcx == "CD") {n1$bcy = "LO"}
        bc4 = "p"; surf = surf + 1
        # Add extra nodes as required
        if (n4$bcx == "CD") {
          nnum = nrow(dfnode) + 1
          n4$nnum = nnum
          n4$bcy = "LO"
          n4 = flags(n4, i+1, j)
          n4 = AirfoilTrans(n4)
          dfnode <- rbind(dfnode, data.frame(
            nnum = nnum,
            x = n4$x, y = n4$y, z = 0, bcx = "padd", bcy = "LO") )
          
          if (n1$bcx == "CD") {
            n1$nnum = nnum - 1
          }
        } else if (n4$bcx == "TE") {
          n1$nnum = nrow(dfnode)
        }
        dfsurf = rbind(dfsurf, data.frame(
          i = i, j = j, surf = surf, enum = enum, side = 4, bc = bc4, nnum1 = n4$nnum, nnum2 = n1$nnum))
      }
    }
    
    # Element Result
    dfrow = data.frame(
      enum    = enum,
      i       = c(i, i+1, i+1, i  ),
      j       = c(j, j  , j+1, j+1),
      nnum    = c(n1$nnum, n2$nnum, n3$nnum, n4$nnum),
      ncorner = c("n1", "n2", "n3", "n4"),
      bc      = c(bc1, bc2, bc3, bc4) )
    dfelem <- rbind(dfelem, dfrow)
  }
}

#--- Plot ----
#--- * Airfoil ----
plotaf <- ggplot(airfoil, aes(x, y)) + geom_path()
plotaf + coord_fixed()
plotaf + 
  geom_point() + 
  geom_point(data = mutate(c, x = -x), colour = "red") +
  coord_fixed(xlim = c(-0.55, -0.45), ylim = c(-0.02, 0.08))
#--- * Nodes ----
plotnode <- ggplot(dfnode, aes(x, y)) +
  geom_path(aes(x, y), airfoil) +
  geom_point(aes(x, y), airfoil, shape = "o") +
  geom_point(aes(colour = nnum)) +
  scale_colour_gradientn(colours = spectralpalette(10)) +
  ggtitle(paste(paste0("NACA", input$NACA), "|", 
                input$alpha, "deg", "|", 
                "polyx =", input$x$polyn, "|",
                "polyy =", input$y$polyn, "|",
                "xn =", input$x$n, "|",
                "yn =", input$y$n, "|"))
plotnode +
  coord_fixed(xlim = c(-0.5, 0.5), ylim = c(-0.3, 0.3))
plotnode +
  coord_fixed(xlim = c(-0.51, -0.40), ylim = c(-0.02, 0.08))

#--- * Surfaces ----
# Create a dataframe of (x,y) to ggplot
plotsurf <- data.frame()
for (i in 1:nrow(dfsurf)) {
  dfrow = data.frame(
    bc   = as.character(dfsurf$bc[i]),
    side = dfsurf$side[i],
    node = c(dfsurf$nnum1[i], dfsurf$nnum2[i]),
    x    = c(filter(dfnode, nnum == dfsurf$nnum1[i])$x, filter(dfnode, nnum == dfsurf$nnum2[i])$x),
    y    = c(filter(dfnode, nnum == dfsurf$nnum1[i])$y, filter(dfnode, nnum == dfsurf$nnum2[i])$y),
    stringsAsFactors = FALSE)
  plotsurf <- rbind(plotsurf, dfrow)
}
plotsurf$side <- as.factor(plotsurf$side)
plotsurf <- plotsurf %>% arrange(node)
plotsurf <- plotsurf %>%  mutate(
  y  = ifelse(node > max(matnnum), y - 0.01*(input$NACA == 0), y),
  bc.side = paste(bc, side, sep = "."))


plotbc <- ggplot(plotsurf, aes(x, y)) +
  geom_path(aes(x, y), airfoil) + 
  geom_line(aes(colour = bc.side)) +
  geom_point(aes(shape = bc.side, colour = bc.side))

plotbc + coord_fixed()
plotbc + coord_fixed(xlim = c(-0.5, 0.5), ylim = c(-0.3, 0.3)) 

#--- * Mesh ----
plotelem <- left_join(dfelem, dfnode, by = "nnum")
plotmesh <- ggplot(plotelem, aes(x, y, group = enum)) +
  geom_path(aes(x, y), airfoil) +
  geom_polygon(aes(colour = enum), fill = NA) +
  geom_point(data = mutate(c, x = -x)) +
  geom_point(data = mutate(c, y = -y)) +
  scale_colour_gradientn(colours = spectralpalette(10))

plotmesh + 
  coord_fixed(expand = FALSE)

plotmesh + 
  coord_fixed(xlim = c(-0.51, -0.40), ylim = c(-0.02, 0.08), expand = FALSE)

plotmesh + 
  coord_fixed(xlim = c(0.40, 0.51), ylim = c(-0.08, 0.02), expand = FALSE)

plotmesh +
  coord_fixed(xlim = c(0.49855, 0.49910), ylim = c(-0.03505, -0.03470))

plotmesh +
  coord_fixed(xlim = c(0.49750, 0.50000), ylim = c(-0.03570, -0.03420))

plotmesh +
  coord_cartesian(xlim = c(-0.4990, -0.498), ylim = c(0.031, 0.038))
