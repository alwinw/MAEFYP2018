#============================#
# Plot Boundary Conditions
# Alwin Wang
#----------------------------#

#--- Set Up                                                       ----
# Use rstudioapi to get saved location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
#--- * Scripts                                                    ----
# Source Required Scripts
source("src_library-manager.R")                                 # Call libraries and install missing ones
source("src_helper-functions.R")                                 # Smaller functions used
# Additional scripts here
# ggplot2 setup (consider moving to a separate script)
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 

D_T       = 0.0005
N_STEP    = 4000

MOD_ALPHA_X = paste0(
  "-(( 5*PI/2*sin(5*PI*(t-0.05)))*heav(t-0.05)+",
  "  (-5*PI/2*sin(5*PI*(t-0.05)))*heav(t-0.25)+",
  "  (-5*PI/2*sin(5*PI*(t-0.05)))*heav(t-0.85)+",
  "  ( 5*PI/2*sin(5*PI*(t-0.05)))*heav(t-1.05))")

u = paste0(
  "( 1/2-1/2*cos(5*PI*(t-0.05)   ))*heav(t-0.05)+",
  "( 1/2-1/2*cos(5*PI*(t-0.05)+PI))*heav(t-0.25)+",
  "(-1/2-1/2*cos(5*PI*(t-0.05)+PI))*heav(t-0.85)+",
  "(-1/2-1/2*cos(5*PI*(t-0.05)   ))*heav(t-1.05)")
x = paste0(
  "-((- 1/40+1/2*t-1/(2*5*PI)*sin(5*PI*(t-0.05)   ))*heav(t-0.05)+",
  "  (- 5/40+1/2*t-1/(2*5*PI)*sin(5*PI*(t-0.05)+PI))*heav(t-0.25)+",
  "  ( 17/40-1/2*t-1/(2*5*PI)*sin(5*PI*(t-0.05)+PI))*heav(t-0.85)+",
  "  ( 21/40-1/2*t-1/(2*5*PI)*sin(5*PI*(t-0.05)   ))*heav(t-1.05))"
)
MOD_ALPHA_X %<>% gsub(MOD_ALPHA_X, "", .) %>%
  gsub("\t", "", .) %>%
  gsub("=", "", .) %>%
  gsub("PI", "pi", .)
eval(parse(text = paste0(
  "BC_", tolower("MOD_ALPHA_X"), " <- function(t) ", MOD_ALPHA_X)),
  envir = .GlobalEnv)
u %<>% gsub(u, "", .) %>%
  gsub("\t", "", .) %>%
  gsub("=", "", .) %>%
  gsub("PI", "pi", .)
eval(parse(text = paste0(
  "BC_", tolower("INFLOW_U"), " <- function(t) ", u)),
  envir = .GlobalEnv)
x %<>% gsub(x, "", .) %>%
  gsub("\t", "", .) %>%
  gsub("=", "", .) %>%
  gsub("PI", "pi", .)
eval(parse(text = paste0(
  "BC_", tolower("AIRFOIL_X"), " <- function(t) ", x)),
  envir = .GlobalEnv)

t = seq(0, D_T*N_STEP, length.out = N_STEP + 1)

plot <- data.frame(
  t = t,
  x = BC_airfoil_x(t),
  u = BC_inflow_u(t),
  a = BC_mod_alpha_x(t)
)

ggplot(plot, aes(t, a)) +
  geom_path() +
  scale_x_continuous(breaks = c(0, seq(0.05, 2.05, 0.2)),
                     minor_breaks = c(0, seq(0.05, 2.05, 0.1)),
                     labels = function(x) ifelse(x==0, "", sprintf("%.2f", x)),
                     limits = c(0, 2)) +
  xlab("Time (s)") +
  ylab(expression("Airfoil Acceleration, a (L / T"^2*")")) +
  scale_y_reverse()
ggsave("BC_airfoil_accel.png", width = 15, height = 6, units = "cm")

ggplot(plot, aes(t, u)) +
  geom_path() +
  scale_x_continuous(breaks = c(0, seq(0.05, 2.05, 0.2)),
                     minor_breaks = c(0, seq(0.05, 2.05, 0.1)),
                     labels = function(x) ifelse(x==0, "", sprintf("%.2f", x)),
                     limits = c(0, 2)) +
  xlab("Time (s)") +
  ylab("Inflow Velocity, V(t) (L / T)")
ggsave("BC_inflow_vel.png", width = 15, height = 6, units = "cm")

ggplot(plot, aes(t, x)) +
  geom_path() +
  scale_x_continuous(breaks = c(0, seq(0.05, 2.05, 0.2)),
                     minor_breaks = c(0, seq(0.05, 2.05, 0.1)),
                     labels = function(x) ifelse(x==0, "", sprintf("%.2f", x)),
                     limits = c(0, 2)) +
  xlab("Time (s)") +
  ylab("Airfoil Displacement, x (L)") +
  scale_y_reverse()
ggsave("BC_airfoil_disp.png", width = 15, height = 6, units = "cm")