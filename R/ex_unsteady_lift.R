#============================#
# Plot Boundary Conditions
# Alwin Wang
#----------------------------#

source("plot_BCs.R")

# t = seq(0, D_T*N_STEP, length.out = N_STEP + 1)
t = seq(0, D_T*N_STEP, length.out = 10)
flow <- data.frame(
  t = t,
  Vinf = BC_inflow_u(t)
)
theta = seq(0, pi, length.out = 10)
alpha = 4*pi/180



# Needs to be a discretised problem ugh. Back to Gamma and x


flx <- read.table(file = "../src-example/NACA0012_remesh_afmc/results/NACA0012_remesh.flx")
colnames(flx) <- c("step", "t", 
                   paste0(c("Fpre", "Fvis", "Ftot"), "x"),
                   paste0(c("Fpre", "Fvis", "Ftot"), "y"),
                   paste0(c("Fpre", "Fvis", "Ftot"), "z"))
flx <- mutate(
  Vinf <- BC_inflow_u()
)

flx$Vinf <- BC_inflow_u(flx$t)  
flx$Gamma <- pi*alpha*flx$Vinf
flx$Gammax <- Gamma/4
