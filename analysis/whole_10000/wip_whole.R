#============================#
# Analysis of Whole Time Domain
# Alwin Wang
#----------------------------#

#--- Set Up                                                       ----
# Use rstudioapi to get saved location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
#--- * Scripts                                                    ----
# Source Required Scripts
srcpath = "../../R/"
source(paste0(srcpath, "src_library-manager.R"))                # Call libraries and install missing ones
source(paste0(srcpath, "src_batch-functions.R"))                # Batch functions used
source(paste0(srcpath, "src_helper-functions.R"))               # Smaller functions used
# Additional scripts here
# ggplot2 setup (consider moving to a separate script)
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 

#--- * Required Paths                                             ----
if (FALSE) {
  #--- > Single File Analysis                                       ----
  saveplot   = "images/"
  airfoil    = "NACA0012r"
  folderpath = "results/"
  seshpath   = "NACA0012r"
  dumppath   = "NACA0012r-249.dump"
  #--- > Required input dataframes                                  ----
  data_airfoil <- data.frame(
    airfoil  = airfoil,
    seshname = seshpath,
    folder   = folderpath,
    seshpath = paste0(folderpath, seshpath),
    stringsAsFactors = FALSE)
  data_mesh <- data.frame(
    data_airfoil,
    tokenword  = "N_P",
    tokenvalue = 8,
    ID         = "NACA0012r-N_P8",
    stringsAsFactors = FALSE)
  data_dump <- data.frame(
    data_mesh,
    dumpfile = dumppath,
    dumppath = paste0(folderpath, dumppath),
    stringsAsFactors = FALSE)
  rm(airfoil, folderpath, seshpath, dumppath)
  #--- > Call Batch Functions                                       ----
  plot = 0
  # Airfoil
  outp_airfoil        <- list(BatchLoadAirfoil(data_airfoil, srcpath))
  names(outp_airfoil) <- data_airfoil$airfoil
  # Mesh
  outp_mesh           <- list(BatchLoadMesh(data_mesh, outp_airfoil, srcpath))
  names(outp_mesh)    <- data_mesh$ID
  # Dump
  outp_dump           <- list(BatchLoadDump(data_dump, outp_mesh, 
                                            plot = c("none"), 
                                            outp = "wall", srcpath))
  names(outp_dump)    <- data_dump$dumppath
  
} else {
  # Batch Test Case
  # Batch list to process
  batchfolder = "results"
  df_batch <- ListSesh(batchfolder, c("airfoil"))
  # List of unique airfoils & N_P
  li_airfoil   <- split (df_batch, df_batch$airfoil)
  li_airfoil   <- lapply(li_airfoil, function(df) df[1,])
  outp_airfoil <- pblapply(li_airfoil, BatchLoadAirfoil, srcpath)
  # List of meshs (should turn this into a function!)
  tokenwords = list("N_P")                                        # Find unique bndry and knot N values
  df_mesh <- df_batch %>%
    rowwise() %>%
    do(data.frame(
      ., LoadSeshTokenWords(.$seshpath, tokenwords), stringsAsFactors = FALSE)) %>%
    filter(!is.na(tokenvalue)) %>% 
    mutate(                                                       # May need to adjust later
      ID = paste0(airfoil))      
  li_mesh   <- split (df_mesh, df_mesh$ID)                        # Unique airfoil and N_P
  li_mesh   <- lapply(li_mesh, function(x) x[1,])                 # Take only the first row in each
  outp_mesh <- pblapply(li_mesh, BatchLoadMesh, outp_airfoil, srcpath)
  # List of dumps
  df_dump <- df_mesh %>% 
    rowwise() %>% 
    do(data.frame(., dumpfile = ListDump(.$folder, .$seshname),
                  stringsAsFactors = FALSE)) %>%
    mutate(dumppath = paste0(folder, dumpfile))
  li_dump   <- split(df_dump, df_dump$dumppath)
  # Use cluster
  cl <- makeCluster(detectCores() - 1)
  # clusterExport(cl, c("outp_mesh", "srcpath"))
  outp_dump <- pblapply(li_dump, BatchLoadDump, outp_mesh, 
                        plot = c("none"), 
                        outp = "wall", srcpath,
                        cl = cl)
  stopCluster(cl)
}

save(outp_dump, file="outp_dump.RData")

#--- > Batch Calc Output                                        ----
outp <- list()
# Wall
outp$wall <- lapply(outp_dump, function(dump) cbind(
  dump$wall, 
  select(dump$data_plot, airfoil, seshname, tokenvalue, ID, time, kinvis, a) ))
outp$wall <- bind_rows(outp$wall)
# BVFa
outp$BVFa <- lapply(outp_dump, function(dump) cbind(
  dump$bvfa,
  select(dump$data_plot, airfoil, seshname, time, kinvis, a) ))
outp$BVFa <- bind_rows(outp$BVFa)
# Circulation
outp$circ <- lapply(outp_dump, function(dump) cbind(
  rbind(select(dump$inte, regn, o),
        data.frame(regn = "W", o = sum(dump$inte[c(1,3,4),2]))),
  select(dump$data_plot, airfoil, seshname, time, kinvis, a) ))
outp$circ <- bind_rows(outp$circ)
# Vortex Impulse
outp$impl <- lapply(outp_dump, function(dump) cbind(
  select(dump$inte, regn, ox, oy),
  select(dump$data_plot, airfoil, seshname, time, kinvis, a) ))
outp$impl <- bind_rows(outp$impl)
outp$impl <- outp$impl %>% 
  arrange(time)
# Vortex Force
outp$forc <- outp$impl %>% 
  filter(regn == "Total") %>% 
  mutate(doxdt = (lead(ox)-ox)/(lead(time)-time),
         doydt = (lead(oy)-oy)/(lead(time)-time) )
# Lift and Drag Force
outp$LD <- read.table(
  paste0(df_batch$seshpath, ".flx"), skip = 3 )
colnames(outp$LD) <- c(
  "step", "time", 
  "Fpre.x", "Fvis.x", "Ftot.x",
  "Fpre.y", "Fvis.y", "Ftot.y",
  "Fpre.z", "Fvis.z", "Ftot.z" )

#--- > Batch Plot Output                                        ----
# BFVa
plot_BVFa <- ggplot(outp$BVFa, aes(time)) +
  annotate("rect", xmin=0.05, xmax=0.25, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  annotate("rect", xmin=0.85, xmax=1.05, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  geom_line(aes(y = RHS, colour = "RHS")) +
  geom_line(aes(y = LHS, colour = "LHS")) +
  xlab("Time (s)") + ylab("Integral of BVF on Airfoil") +
  scale_x_continuous(breaks = c(0, seq(0.05, 2.5, 0.2), 2.5),
                     minor_breaks = c(0, seq(0.05, 2.5, 0.1), 2.5),
                     labels = function(x) ifelse(x==0 | x==2.5, "", sprintf("%.2f", x)),
                     limits = c(0, 2.5))
ggsave(paste0("BVFa.png"), plot_BVFa,
       scale = 2, width = 10, height = 4, units = "cm", dpi = 300)
# Circulation
plot_circ <- ggplot(outp$circ %>% filter(regn != "W"), 
       aes(time, o, group = regn, colour = regn, linetype = regn)) +
  annotate("rect", xmin=0.05, xmax=0.25, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  annotate("rect", xmin=0.85, xmax=1.05, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  geom_line() +
  xlab("Time (s)") + ylab("Circulation") +
  scale_x_continuous(breaks = c(0, seq(0.05, 2.5, 0.2), 2.5),
                     minor_breaks = c(0, seq(0.05, 2.5, 0.1), 2.5),
                     labels = function(x) ifelse(x==0 | x==2.5, "", sprintf("%.2f", x)),
                     limits = c(0, 2.5))
ggsave(paste0("circ.png"), plot_circ,
       scale = 2, width = 10, height = 6, units = "cm", dpi = 300)
# Vortex Impulse
plot_impl <- ggplot(outp$impl %>% filter(regn == "Total"), 
       aes(time, group = regn)) +
  annotate("rect", xmin=0.05, xmax=0.25, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  annotate("rect", xmin=0.85, xmax=1.05, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  geom_line(aes(y = ox, colour = "ox")) +
  geom_line(aes(y = oy, colour = "oy")) +
  xlab("Time (s)") + ylab("Vortex Impulse") +
  scale_x_continuous(breaks = c(0, seq(0.05, 2.5, 0.2), 2.5),
                     minor_breaks = c(0, seq(0.05, 2.5, 0.1), 2.5),
                     labels = function(x) ifelse(x==0 | x==2.5, "", sprintf("%.2f", x)),
                     limits = c(0, 2.5))
ggsave(paste0("impl.png"), plot_impl,
       scale = 2, width = 10, height = 6, units = "cm", dpi = 300)
# Vortex Force
plot_forc <- ggplot(outp$forc, 
       aes(time, group = regn)) +
  annotate("rect", xmin=0.05, xmax=0.25, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  annotate("rect", xmin=0.85, xmax=1.05, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  geom_line(aes(y = doxdt, colour = "doxdt")) +
  geom_line(aes(y = doydt, colour = "doydt")) +
  xlab("Time (s)") + ylab("Vortex Impulse") +
  scale_x_continuous(breaks = c(0, seq(0.05, 2.5, 0.2), 2.5),
                     minor_breaks = c(0, seq(0.05, 2.5, 0.1), 2.5),
                     labels = function(x) ifelse(x==0 | x==2.5, "", sprintf("%.2f", x)),
                     limits = c(0, 2.5))
ggsave(paste0("forc.png"), plot_forc,
       scale = 2, width = 10, height = 6, units = "cm", dpi = 300)
# Lift and Drag
plot_LD <- ggplot(outp$LD, aes(time)) +
  annotate("rect", xmin=0.05, xmax=0.25, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  annotate("rect", xmin=0.85, xmax=1.05, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  geom_line(aes(y = Ftot.x, colour = "Ftot.x")) +
  geom_line(aes(y = Ftot.y, colour = "Ftot.y"))
ggsave(paste0("LD.png"), plot_LD,
       scale = 2, width = 10, height = 6, units = "cm", dpi = 300)

# Plot N-S
plot_t = 0.05


PlotNS <- function(plot_t) {
  plot_ns    <- filter(outp$wall, time == plot_t)
  
  plot_setup <- PlotSetup(plot_ns, plot_ns[1,])
  plot_setup$title <- paste0(
    paste("Time:",                sprintf("%05.3f",  plot_ns$time  )), "   ",
    paste("Acceleration:",        sprintf("%+07.4f", plot_ns$a   )) )
  
  plot_nstheme <- ggplot(plot_ns, aes(s)) + 
    geom_vline(xintercept = as.numeric(plot_setup$vlines),        # Vertical lines for LE, TE, LE
               colour = "grey", linetype = "dashed") +
    geom_label(aes(x, y, label = labels), plot_setup$surf,        # Surface labels
               # geom_label(aes(x, y+20, label = labels), plot_setup$surf,        # Surface labels
               colour = "grey") +
    xlab("s") + 
    scale_x_continuous(breaks = plot_setup$xbreaks, 
                       labels = function(x) sprintf("%.2f", x)) +
    ylab(NULL) + ylim(c(-40, 30)) +
    # ylab(NULL) + ylim(c(-20, 20)) +
    scale_color_manual(
      name = "Legend",
      values = c("dp/ds" = "red", "dV/dt" = "blue", "LHS" = "purple", "RHS" = "purple"),
      labels = c(
        expression(-frac(1, rho)~frac(partialdiff*p, partialdiff*s)), 
        expression(-frac(partialdiff*V, partialdiff*t)), 
        expression(-bgroup("(",
                           frac(1, rho)~frac(partialdiff*p, partialdiff*s) +
                             frac(partialdiff*V, partialdiff*t), ")")),
        expression(-nu~frac(partialdiff*omega, partialdiff*z))),
      guide = guide_legend(
        override.aes = list(
          linetype = c(rep("solid", 2), "dashed", "blank"),
          # shape = c(rep(NA, 3), "O"),
          shape = c(rep(NA, 3), 20),
          alpha = rep(1, 4)))) +
    theme(legend.key.size = unit(2.25, "lines"),
          legend.text.align = 0.5,
          legend.direction = "vertical", 
          legend.position = "right",
          legend.background = element_rect(colour = "black", size = 0.3),
          plot.title = element_text(size = 10))
  
  plot_nsG <- plot_nstheme +
    geom_path(aes(y = -as, colour = "dV/dt")) +                   # Acceleration terms
    geom_path(aes(y = + dpdsG, colour = "dp/ds")) +               # Pressure field
    geom_path(aes(y = - as + dpdsG, colour = "LHS"),              # LHS, acceleration + pressure
              linetype = "dashed") +
    geom_point(aes(s, RHSG, colour = "RHS"), shape = 'o', alpha = 0.3) +
    # geom_point(aes(s, -dodzG*plot_data$kinvis, colour = "RHS"),    # RHS, v * dw/dz
               # plot_offs, shape = 20) +
    # plot_offs, shape = "o", alpha = 0.3) +
    ggtitle(plot_setup$title)
  
  print(plot_nsG)
  
  save = paste0("Rslt_Re", 
                sprintf("%05d", round(1/plot_ns$kinvis[1])), 
                "-NS-t", 
                sprintf("%0.2f", plot_t),
                ".png")
  
  ggsave(save, plot_nsG,
         scale = 2, width = 10, height = 4.5, units = "cm", dpi = 200)
  return(save)
}

for (i in c(0.05, 0.10, 0.15, 0.20,
            0.25, 0.45, 0.65, 0.85,
            0.90, 0.95, 1.00, 1.05,
            1.25, 1.45, 1.65, 1.85)) {
  PlotNS(i)
}
