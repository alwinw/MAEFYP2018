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
  polyn = 3
)

#--- Polynomial ----
c = data.frame(x = 0.5*cos(input$alpha*pi/180),
               y = 0.5*sin(input$alpha*pi/180))
# x
dx = (input$x$max - input$x$min)/input$x$n
listx = list(
  LE = seq(input$x$min, -c$x, length.out = round((-c$x - input$x$min)/dx) + 1),
  CD = seq(-c$x, c$x, length.out = 3*round(c$x*2/dx) + 1),
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
# Node numbering
nx = nrow(dfx)
ny = nrow(dfy)
matnnum <- matrix(1:(nx*ny), nrow = ny, ncol = nx, byrow = TRUE)
dfnode <- data.frame()
for (j in 1:ny) {
  for (i in 1:nx) {
    dfrow = data.frame(
      nnum = matnnum[j, i], x = dfx$x[i], y = dfy$y[j], z = 0, bcx = dfx$bcx[i])
    dfrow = flags(dfrow, i, j)
    dfnode = rbind(dfnode, dfrow)
  } 
}

ggplot(dfnode, aes(x, y, colour = nnum)) +
  geom_point() + 
  scale_colour_gradientn(colours = spectralpalette(10))

# Elements
dfelem <- data.frame()
enum = 0
dfsurf <- data.frame()
surf = 0
for (j in 1:(ny - 1)) {
  for (i in 1:(nx - 1)) {
    # 1 ---- 4
    # | enum |
    # 2 ---- 3
    enum = enum + 1
    # Determine node corners
    n1 = data.frame(x = dfx$x[i  ], y = dfy$y[j  ], nnum = matnnum[j  , i  ], bcx = dfx$bcx[i  ], bcy = dfy$bcy[j  ])
    n2 = data.frame(x = dfx$x[i  ], y = dfy$y[j+1], nnum = matnnum[j+1, i  ], bcx = dfx$bcx[i  ], bcy = dfy$bcy[j+1])
    n3 = data.frame(x = dfx$x[i+1], y = dfy$y[j+1], nnum = matnnum[j+1, i+1], bcx = dfx$bcx[i+1], bcy = dfy$bcy[j+1])
    n4 = data.frame(x = dfx$x[i+1], y = dfy$y[j  ], nnum = matnnum[j  , i+1], bcx = dfx$bcx[i+1], bcy = dfy$bcy[j  ])
    # Determine y values based on bcx
    n1 = flags(n1, i  , j  )
    n2 = flags(n2, i  , j+1)
    n3 = flags(n3, i+1, j+1)
    n4 = flags(n4, i+1, j  )
    # Determine wall bcs
    bc1 = NA; bc2 = NA; bc3 = NA; bc4 = NA;
    if (j == 1) { # TOP | side 4
      bc4 = "s"; surf = surf + 1
      dfsurf = rbind(dfsurf, data.frame(
        i = i, j = j, surf = surf, enum = enum, side = 4, bc = bc4, nnum1 = n4$nnum, nnum2 = n1$nnum)) }
    if (i == (nx - 1)) { # RIGHT | side 3
      bc3 = "o"; surf = surf + 1
      dfsurf = rbind(dfsurf, data.frame(
        i = i, j = j, surf = surf, enum = enum, side = 3, bc = bc3, nnum1 = n3$nnum, nnum2 = n4$nnum)) }
    if (j == (ny - 1)) { # BOTTOM | side 2
      bc2 = "w"; surf = surf + 1
      dfsurf = rbind(dfsurf, data.frame(
        i = i, j = j, surf = surf, enum = enum, side = 2, bc = bc2, nnum1 = n2$nnum, nnum2 = n3$nnum)) }
    if (i == 1) { # LEFT | side 1
      bc1 = "v"; surf = surf + 1
      dfsurf = rbind(dfsurf, data.frame(
        i = i, j = j, surf = surf, enum = enum, side = 1, bc = bc1, nnum1 = n1$nnum, nnum2 = n2$nnum)) }

    # Side 3 airfoil bc
    if (n2$bcy == "CD" & n3$bcy == "CD") { # Element above plate | side 2
      if (n3$bcx %in% c("CD", "TE")) {
        bc2 = "p"; surf = surf + 1
        dfsurf = rbind(dfsurf, data.frame(
          i = i, j = j, surf = surf, enum = enum, side = 2, bc = bc2, nnum1 = n2$nnum, nnum2 = n3$nnum)) }
    }
    # Side 1 airfoil bc
    if (n4$bcy == "CD" & n1$bcy == "CD") {
      if (n4$bcx %in% c("CD", "TE")) { # Element below plote | side 4
        bc4 = "p"; surf = surf + 1
        # Add extra nodes as required
        if (n4$bcx == "CD") {
          nnum = nrow(dfnode) + 1
          dfnode <- rbind(dfnode, data.frame(
            nnum = nnum,
            x = n4$x, y = n4$y, z = 0, bcx = "padd") )
          n4$nnum = nnum
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
    
    # Result
    dfrow = data.frame(
      enum    = enum,
      i       = c(i, i+1, i+1, i  ),
      j       = c(j, j  , j+1, j+1),
      nnum    = c(n1$nnum, n2$nnum, n3$nnum, n4$nnum),
      ncorner = c("n1", "n2", "n3", "n4"),
      x       = c(n1$x, n2$x, n3$x, n4$x),
      y       = c(n1$y, n2$y, n3$y, n4$y),
      bc      = c(bc1, bc2, bc3, bc4) )
    dfelem <- rbind(dfelem, dfrow)
  }
}

plotmesh <- ggplot(dfelem, aes(x, y, group = enum)) +
  geom_polygon(aes(colour = enum), fill = NA) +
  # geom_path(aes(colour = bc, group = bc), dfelem %>% filter(!is.na(bc))) +
  geom_point(data = mutate(c, x = -x)) +
  geom_point(data = mutate(c, y = -y)) + 
  scale_colour_gradientn(colours = spectralpalette(10))
  
plotmesh + coord_fixed()
plotmesh + coord_cartesian(xlim = c(-0.6,  0.6), ylim = c(-0.1, 0.1))
plotmesh + coord_cartesian(xlim = c(-0.6, -0.4), ylim = c(-0.1, 0.1))
plotmesh + coord_cartesian(xlim = c( 0.4,  0.6), ylim = c(-0.1, 0.1))

# Double check the BCs are correct
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
  y  = ifelse(node > max(matnnum), y - 0.01, y))

plotsurf <- rbind(plotsurf %>% filter(node <= max(matnnum)),
                  plotsurf %>% filter(node > max(matnnum)) %>% arrange(-node))


plotbc <- ggplot(plotsurf, aes(x, y, group = bc, colour = bc)) + 
  geom_polygon(fill = NA) +
  geom_point(aes(shape = bc))

plotbc + coord_fixed()
plotbc + coord_fixed(xlim = c(-0.6,  0.6), ylim = c(-0.1, 0.1))

plotside <- ggplot(plotsurf, aes(x, y, group = bc, colour = side)) + 
  geom_polygon(fill = NA) +
  geom_point(aes(shape = side))

plotside + coord_fixed()
plotside + coord_fixed(xlim = c(-0.6,  0.6), ylim = c(-0.1, 0.1))


#--- Print the Output ----
name = "plate.sesh"
# Header
cat(paste0(
  "<USER>\n", 
  "  u = 0.0\n",
  "  v = 0.0\n",
  "  p = 0.0\n",
  "</USER>\n",
  "\n",
  "<FIELDS>\n",
  "  u v p\n",
  "</FIELDS>\n",
  "\n",
  "<TOKENS>\n",
  "  KINVIS    = 1./10000.\n",
  "  \n",
  "  D_T       = 0.0001 \n",
  "  N_STEP    = 2000\n",
  "  N_TIME    = 2\n",
  "  \n",
  "  N_P       = 8\n",
  "  N_Z       = 1\n",
  "  LZ        = 1.0\n",
  "  BETA      = TWOPI/LZ\n",
  "  \n",
  "  IO_CFL    = 50\n",
  "  IO_FLD    = 10\n",
  "  IO_HIS    = 1\n",
  "  \n",
  "  AVERAGE   = 0\n",
  "  CHKPOINT  = 0\n",
  "  ITERATIVE = 1\n",
  "  TBCS      = 1\n",
  "</TOKENS>\n",
  "\n",
  "<FORCE>\n",
  "  MOD_A_X = 1\n",
  "  MOD_A_Y = 0\n",
  "  MOD_ALPHA_X = -((5*PI/2*sin(5*PI*(t-0.05)))*heav(t-0.05)+(-5*PI/2*sin(5*PI*(t-0.05)))*heav(t-0.25)+(-5*PI/2*sin(5*PI*(t-0.05)))*heav(t-0.85)+(5*PI/2*sin(5*PI*(t-0.05)))*heav(t-1.05)) \n",
  "  MOD_ALPHA_Y = 0\n",
  "</FORCE>\n",
  "  \n",
  "<GROUPS NUMBER=5>\n",
  "  1     o     outflow\n",
  "  2     p     wall\n",
  "  3     w     group_name\n",
  "  4     s     group_name\n",
  "  5     v     inflow\n",
  "</GROUPS>\n",
  "  \n",
  "<BCS NUMBER=5>\n",
  "  1     o     3\n",
  "              <N>     u = 0      </N>\n",
  "              <N>     v = 0      </N>\n",
  "              <D>     p = 0      </D>\n",
  "  2     p     3\n",
  "              <D>     u = 0      </D>\n",
  "              <D>     v = 0      </D>\n",
  "              <H>     p = 0      </H>\n",
  "  3     w     3\n",
  "              <N>     u = 0      </N>\n",
  "              <D>     v = 0      </D>\n",
  "              <H>     p = 0      </H>\n",
  "  4     s     3\n",
  "              <N>     u = 0      </N>\n",
  "              <D>     v = 0      </D>\n",
  "              <H>     p = 0      </H>\n",
  "  5     v     3\n",
  "              <D>     u = (1/2-1/2*cos(5*PI*(t-0.05)))*heav(t-0.05)+(1/2-1/2*cos(5*PI*(t-0.05)+PI))*heav(t-0.25)+(-1/2-1/2*cos(5*PI*(t-0.05)+PI))*heav(t-0.85)+(-1/2-1/2*cos(5*PI*(t-0.05)))*heav(t-1.05)   </D>\n",
  "              <D>     v = 0      </D>\n",
  "              <H>     p = 0      </H>\n",
  "</BCS>\n",
  "\n"),
  file = name, append = FALSE)
# Nodes
cat(paste0("<NODES NUMBER=", nrow(dfnode), ">\n"), file = name, append = TRUE)
for (n in 1:nrow(dfnode)) {
  cat(sprintf("%6d%16.10f%16.10f%16.10f\n", dfnode[n, 1], dfnode[n,2], dfnode[n,3], dfnode[n,4]), 
      file = name, append = TRUE)
}
cat("</NODES>\n\n", file = name, append = TRUE)
# Elements
cat(paste0("<ELEMENTS NUMBER=", nrow(dfelem)/4, ">\n"), file = name, append = TRUE)
for (n in 1:(nrow(dfelem)/4)) {
  cat(sprintf("%6d   <Q>%6d%6d%6d%6d   </Q>\n", dfelem[4*n, 1], 
              dfelem[4*n-3,4], dfelem[4*n-2,4], dfelem[4*n-1,4], dfelem[4*n,4]),
      file = name, append = TRUE)
}
cat("</ELEMENTS>\n\n", file = name, append = TRUE)
# Surfaces
dfsurf <- dfsurf %>% arrange(bc, enum, side)
cat(paste0("<SURFACES NUMBER=", nrow(dfsurf), ">\n"), file = name, append = TRUE)
for (n in 1:nrow(dfsurf)) {
  cat(sprintf("%6d%6d%6d   <B>%9s      </B>\n", dfsurf[n, 3], dfsurf[n,4], dfsurf[n,5], dfsurf[n,6]), 
      file = name, append = TRUE)
}
cat("</SURFACES>", file = name, append = TRUE)

