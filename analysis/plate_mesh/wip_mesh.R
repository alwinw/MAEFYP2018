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
  y = data.frame(min = -1.25, max = 1.25,   n = 51),
  alpha = 4, # degrees
  polyn = 2
)

#--- Calculation ----
c = data.frame(x = 0.5*cos(input$alpha*pi/180),
               y = 0.5*sin(input$alpha*pi/180))
# x
dx = (input$x$max - input$x$min)/input$x$n
listx = list(
  LE = seq(input$x$min, -c$x, length.out = round((-c$x - input$x$min)/dx) + 1),
  CD = seq(-c$x, c$x, length.out = 1.5*round(c$x*2/dx) + 1),
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
ggplot() +
  # geom_vline(xintercept = listx$LE, linetype = "dotted") +
  # geom_vline(xintercept = listx$CD, linetype = "dotted") +
  # geom_vline(xintercept = listx$TE, linetype = "dotted") +
  geom_vline(xintercept = outpx$LE, colour = "blue") +
  geom_vline(xintercept = outpx$CD, colour = "red") +
  geom_vline(xintercept = outpx$TE, colour = "blue") +
  geom_vline(xintercept = outpx$WK, colour = "green")

# y
dy = (input$y$max - input$y$min)/input$y$n
listy = list(
  LO = seq(input$y$min, -c$y, length.out = round((-c$y - input$y$min)/dy) + 1),
  CD = seq(-c$y, c$y, length.out = 1.5*round(c$y*2/dy) + 1),
  UP = seq(c$y, -input$y$min, length.out = round((-input$y$min - c$y)/dy) + 1) )

polyy = list(
  LO = polypow(c(1,  c$x), input$polyn + 1), 
  CD = polyint(polymul(polypow(c(1,-c$x),input$polyn), polypow(c(1,c$x),input$polyn))),
  UP = polypow(c(1, -c$x), input$polyn + 1) )
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

ggplot() +
  geom_hline(yintercept = listy$LO, linetype = "dotted") +
  geom_hline(yintercept = listy$CD, linetype = "dotted") +
  geom_hline(yintercept = listy$UP, linetype = "dotted")
