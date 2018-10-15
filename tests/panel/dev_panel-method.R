#============================#
# Panel Method
# Alwin Wang
#----------------------------#

# Panel method: Vortex Panels of Linearly Varying Strength

#--- Set Up                                                       ----
# Use rstudioapi to get saved location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
#--- * Scripts                                                    ----
# Source Required Scripts
srcpath = "../../R/"
source(paste0(srcpath, "src_library-manager.R"))                # Call libraries and install missing ones
source(paste0("../../analysis/panel/", "src_panel-functions.R"))
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 


#--- Test Cases                                                   ----
#--- * Diamond Airfoil                                            ----
# Source: MAE4409 Lecture notes; Hugh Blackburn's C code
coord <- data.frame(
  x = c(1,  0.3,  0,  0.3, 1),
  y = c(0, -0.05, 0, 0.05, 0))
panels <- VLPanel(coord)
ggplot(panels) +
  geom_segment(aes(x1, y1, xend = x, yend = y),
               arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
  geom_segment(aes(x, y, xend = x2, yend = y2)) +
  geom_segment(aes(x, y, xend = x+nx, yend = y+ny)) +
  geom_point(aes(x, y), coord) +
  geom_point(aes(x, y), coord[1,], colour = "red") +
  coord_fixed()
VLSteady(coord, aoa = 5, U = 1)
actual <- data.frame(gamma = c(-0.8862, -0.9905, 0.4848, 1.233, 0.8862))

#--- * NACA2412 Airfoil                                           ----
# Source: Kuethe & Chow pg 121
coord <- data.frame(
  x = c(1.000, 0.933, 0.750, 0.500, 0.250, 0.067, 0.000, 0.067, 0.250, 0.500, 0.750, 0.933, 1.000),
  y = c(0.000,-0.005,-0.017,-0.033,-0.042,-0.033, 0.000, 0.045, 0.076, 0.072, 0.044, 0.013, 0.000) )
panels <- VLPanel(coord)
ggplot(coord, aes(x, y)) + geom_path() + geom_point() + 
  coord_fixed() +
  geom_segment(aes(x, y, xend = x+nx, yend = y+ny), panels)
VLSteady(coord, aoa = 8, U = 1)

#--- * NACA0012 Airfoil                                           ----
# Source: Hugh Blackburn
coord <- data.frame(
  x = c(1.000000, 0.904508, 0.654508, 0.345491, 0.095491,
        0.000000, 0.095491, 0.345491, 0.654508, 0.904508,
        1.000000),
  y = c(0.000000,-0.013914,-0.040917,-0.059575,-0.046049,
        0.000000, 0.046049, 0.059575, 0.040917, 0.013914,
        0.000000) )
panels <- VLPanel(coord)
ggplot(coord, aes(x, y)) + geom_path() + geom_point() + 
  coord_fixed() +
  geom_segment(aes(x, y, xend = x+nx, yend = y+ny), panels)
VLSteady(coord, aoa = 4, U = 1)

#--- * Cylinder                                                   ----
# Source: Katz & Plotkin pg 65
R = 3; Uinf = 20
t = seq(0, -2*pi, length.out = 201)
coord <- data.frame(
  t = t, x = R*cos(t), y = R*sin(t) )
panels <- VLPanel(coord)
ggplot(panels) +
  geom_segment(aes(x1, y1, xend = x, yend = y),
               arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
  geom_segment(aes(x, y, xend = x2, yend = y2)) +
  geom_segment(aes(x, y, xend = x+nx, yend = y+ny)) +
  geom_point(aes(x, y), coord) +
  geom_point(aes(x, y), coord[1,], colour = "red") +
  coord_fixed()
steady <- VLSteady(coord, aoa = 0, U = Uinf)
actual <- data.frame(t = seq(-pi, pi, length.out = 101)) %>% 
  mutate(
    cp  = 1 - 4*sin(t)^2,
    phi = Uinf*cos(t)*2*R )
# Coefficient of Pressure
ggplot() +
  geom_point(aes(atan2(y,x), cp), 
            steady$i, colour = "red") +
  geom_line(aes(ifelse(t < -pi, t+2*pi, t), cp), 
            actual, colour = "blue") +
  scale_y_reverse() +
  ggtitle("Coefficient of Pressure on a Cylinder")
# Velocity potential
ggplot() +
  geom_point(aes(atan2(y,x), phi), 
             steady$i, colour = "red") +
  geom_line(aes(ifelse(t < -pi, t+2*pi, t), phi), 
            actual, colour = "blue") +
  ggtitle("Velocity Potential")
# Determine velocity potential manually by integration
ggplot() +
  geom_point(aes(atan2(y,x), gamma), 
             steady$j, colour = "red") +
  ggtitle("Vorticity Distribution")

#--- * Joukowski Airfoil                                          ----
# Theory of Wing Sections
# Theoretical and applied aerodynamics Ch2 
