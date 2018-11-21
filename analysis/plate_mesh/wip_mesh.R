#============================#
# Mesh Generation
# Alwin Wang
#----------------------------#

#--- Set Up                                                       ----
# Use rstudioapi to get saved location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
#--- * Scripts                                                    ----
# Source Required Scripts
srcpath = "../../R/"
source(paste0(srcpath, "src_library-manager.R"))                # Call libraries and install missing ones
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 

#--- Inputs ----
input <- list(
  x = data.frame(min = -1.75, max = 5.4375, n = 51),
  y = data.frame(min = -1.25, max = 1.25,   n = 21),
  alpha = 4, # degrees
  polyn = 2
)

#--- Polynomial ----
c = data.frame(x = 0.5*cos(input$alpha*pi/180),
               y = 0.5*sin(input$alpha*pi/180))
# x
dx = (input$x$max - input$x$min)/input$x$n
listx = list(
  LE = seq(input$x$min, -c$x, length.out = round((-c$x - input$x$min)/dx) + 1),
  CD = seq(-c$x, c$x, length.out = 2*round(c$x*2/dx) + 1),
  TE = seq(c$x, -input$x$min, length.out = round((-input$x$min - c$x)/dx) + 1) )
polyx = list(
  LE = polypow(c(1,  c$x), input$polyn + 1), 
  CD = polyint(polymul(polypow(c(1,-c$x),input$polyn), polypow(c(1,c$x),input$polyn))),
  TE = polypow(c(1, -c$x), input$polyn + 1) )
polyx$LE = polyadd(polyx$LE/((input$x$min + c$x)^input$polyn), c(-c$x))
polyx$CD = c$x/(polyval(polyx$CD, c$x)) * polyx$CD
polyx$TE = polyadd(polyx$TE/((-input$x$min - c$x)^input$polyn), c( c$x))
outpx = list(
  LE = polyval(polyx$LE, listx$LE),
  CD = polyval(polyx$CD, listx$CD),
  TE = polyval(polyx$TE, listx$TE) )
outpx$WK = seq(outpx$TE[length(outpx$TE)], input$x$max, 
               length.out = round((input$x$max - outpx$TE[length(outpx$TE)])/
                                    (outpx$TE[length(outpx$TE)] - outpx$TE[length(outpx$TE)-1])))

# y
listy = list(
  LO = seq(input$y$min, c$y, length.out = round(input$y$n/2) + 1),
  UP = seq(c$y, input$y$max, length.out = round(input$y$n/2) + 1) )

polyy = list(
  LO = polypow(c(1, -c$y), input$polyn + 1), 
  UP = polypow(c(1, -c$y), input$polyn + 1) )
polyy$LO = polyadd(polyy$LO/((input$y$max + c$y)^input$polyn)*(-2*mod(input$polyn, 2) + 1), c(c$y))
polyy$UP = polyadd(polyy$UP/((input$y$max - c$y)^input$polyn), c(c$y))

outpy = list(
  LO = polyval(polyy$LO, listy$LO),
  UP = polyval(polyy$UP, listy$UP) )

# Plots
# ggplot() +
#   geom_vline(xintercept = listx$LE, linetype = "dotted") +
#   geom_vline(xintercept = listx$CD, linetype = "dotted") +
#   geom_vline(xintercept = listx$TE, linetype = "dotted") +
#   geom_hline(yintercept = listy$LO, linetype = "dotted") +
#   geom_hline(yintercept = listy$UP, linetype = "dotted")

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

#--- Output ----
# Lists
dfx <- data.frame(
  x   = c(outpx$LE, outpx$CD[-1], outpx$TE[-1], outpx$WK[-1]),
  bcx = c("v" , rep("FD", length(outpx$LE) - 2), 
          "LE", rep("CD", length(outpx$CD) - 2),
          "TE", rep("BD", length(outpx$TE) - 2),
                rep("BD", length(outpx$WK) - 1),
          "o") )
dfy <- data.frame(
  y   = c(outpy$LO, outpy$UP[-1]),
  bcy = c("w" , rep("LO", length(outpy$LO) - 2),
          "CD", rep("UP", length(outpy$UP) - 2),
          "s") )
dfy = dfy[nrow(dfy):1,]
dfyd <- data.frame(
  y   = -rev(dfy$y),
  bcy = dfy$bcy )
# Node numbering
nx = nrow(dfx)
ny = nrow(dfy)
matnnum <- matrix(1:(nx*ny), nrow = ny, ncol = nx, byrow = TRUE)
# flags
# FD (use dfy) | BD (use rev(dfy))
# CD (need to interpolate between dfy and rev(dfy))
# UP (use original nodes) | LO (add extra for wall)
flags <- function(n1, i, j) {
  if (n1$bcx == "CD") {
    n1$y = (dfyd$y[j] - n1$y)/(c$x*2)*(n1$x + c$x) + n1$y
  } else if (n1$bcx %in% c("TE", "BD", "o")) {
    n1$y = dfyd$y[j]
  }
  return(n1)
}

bc = function(bcm, bcn) {
  if (bcm == bcn) {
    bc = as.character(bcm)
    if (bc %in% c("FD", "BD", "LO", "UP")) bc = NA
    if (bc %in% c("LE", "TE", "CD")) bc = "p"
  } else {
    bc = NA
  }
  return(bc)
}

# Elements
dfelem <- data.frame()
enum = 0
for (j in 1:(ny - 1)) {
  for (i in 1:(nx - 1)) {
    enum = enum + 1
    # Determine node corners
    n1 = data.frame(x = dfx$x[i  ], y = dfy$y[j  ], nnum = matnnum[j  , i  ], bcx = dfx$bcx[i  ], bcy = dfy$bcy[j  ])
    n2 = data.frame(x = dfx$x[i+1], y = dfy$y[j  ], nnum = matnnum[j  , i+1], bcx = dfx$bcx[i+1], bcy = dfy$bcy[j  ])
    n3 = data.frame(x = dfx$x[i+1], y = dfy$y[j+1], nnum = matnnum[j+1, i+1], bcx = dfx$bcx[i+1], bcy = dfy$bcy[j+1])
    n4 = data.frame(x = dfx$x[i  ], y = dfy$y[j+1], nnum = matnnum[j+1, i  ], bcx = dfx$bcx[i  ], bcy = dfy$bcy[j+1])
    # Determine y values based on bcx
    n1 = flags(n1, i  , j  )
    n2 = flags(n2, i+1, j  )
    n3 = flags(n3, i+1, j+1)
    n4 = flags(n4, i  , j+1)
    # Determine bcs
    bc1 = bc(n1$bcy, n2$bcy)
    bc2 = bc(n2$bcx, n3$bcx)
    bc3 = bc(n3$bcy, n4$bcy)
    bc4 = bc(n4$bcx, n1$bcx)
    
    dfrow = data.frame(
      enum    = enum,
      ncorner = c("n1", "n2", "n3", "n4"),
      x       = c(n1$x, n2$x, n3$x, n4$x),
      y       = c(n1$y, n2$y, n3$y, n4$y),
      bc      = c(bc1, bc2, bc3, bc4)
    )
    dfelem <- rbind(dfelem, dfrow)
  }
}

plotmesh <- ggplot(dfelem, aes(x, y, group = enum)) +
  geom_polygon(aes(), colour = "grey90", fill = NA) +
  # geom_path(aes(colour = bc, group = bc), filter(dfelem, !is.na(bc))) +
  geom_point(data = mutate(c, x = -x)) +
  geom_point(data = mutate(c, y = -y))
  
plotmesh + coord_fixed()
plotmesh + coord_cartesian(xlim = c(-0.6, 0.6), ylim = c(-0.1, 0.1))
plotmesh + coord_cartesian(xlim = c(-0.6, -0.4), ylim = c(-0.1, 0.1))
