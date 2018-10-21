#============================#
# Convergence Analysis
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
  folderpath = "NACA0012r/results/NACA0012r-N_P03/"
  seshpath   = "naca0012r-N_P03"
  dumppath   = "naca0012r-N_P03-00.dump"
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
    tokenvalue = 3,
    ID         = "NACA0012r-N_P3",
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
  outp_dump           <- list(BatchLoadDump(data_dump, outp_mesh, plot, "node", srcpath))
  names(outp_dump)    <- data_dump$dumppath
} else {
  # Batch Test Case
}

#--- * Convergence Analysis                                       ----
#--- > Batch Calculation                                          ----
# Batch list to process
batchfolder = "."
df_batch <- ListSesh(batchfolder, c("airfoil", "a", "b", "seshname"))
df_batch$airfoil = df_batch$seshname                            # Airfoil changes per N_P
# List of unique airfoils & N_P
li_airfoil   <- split (df_batch, df_batch$airfoil)
li_airfoil   <- lapply(li_airfoil, function(df) df[1,])
outp_airfoil <- pblapply(li_airfoil, BatchLoadAirfoil, srcpath)
# List of meshs (should turn this into a function!)
tokenwords = list("N_P")                                        # Find unique bndry and knot N values
df_mesh <- df_batch %>%
  rowwise() %>%
  do(data.frame(., LoadSeshTokenWords(.$seshpath, tokenwords),
                stringsAsFactors = FALSE)) %>%
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
# outp_dump <- pblapply(li_dump, BatchLoadDump, outp_mesh, 0, "node", srcpath)
# Use cluster
cl <- makeCluster(detectCores() - 1)
# clusterExport(cl, c("outp_mesh", "srcpath"))
outp_dump <- pblapply(li_dump, BatchLoadDump, outp_mesh, 
                      plot = c("none"), 
                      outp = "node", srcpath,
                      cl = cl)
stopCluster(cl)
#--- > Batch Calc Output                                        ----
outp <- list()
# Wall
outp$wall <- lapply(outp_dump, function(dump) cbind(dump$wall, 
  select(dump$data_plot, airfoil, seshname, tokenvalue, ID, time, kinvis, a)))
outp$wall <- bind_rows(outp$wall) %>% 
  separate(airfoil, c("airfoilname", "junk"), sep = "-", remove = FALSE) %>% 
  select(-junk)
# Node
outp$node <- lapply(outp_dump, function(dump) cbind(dump$node, 
  select(dump$data_plot, airfoil, seshname, tokenvalue, ID, time, kinvis, a)))
outp$node <- bind_rows(outp$node) %>% 
  separate(airfoil, c("airfoilname", "junk"), sep = "-", remove = FALSE) %>% 
  select(-junk)
# BVFa
outp$BVFa <- lapply(outp_dump, function(dump) cbind(
  dump$bvfa,
  select(dump$data_plot, airfoil, seshname, tokenvalue, ID, time, kinvis, a)))
outp$BVFa <- bind_rows(outp$BVFa) %>% 
  separate(airfoil, c("airfoilname", "junk"), sep = "-", remove = FALSE) %>% 
  select(-junk)
# Inte
outp$inte <- lapply(outp_dump, function(dump) cbind(
  dump$inte,
  select(dump$data_plot, airfoil, seshname, tokenvalue, ID, time, kinvis, a)))
outp$inte <- bind_rows(outp$inte) %>% 
  separate(airfoil, c("airfoilname", "junk"), sep = "-", remove = FALSE) %>% 
  select(-junk)

#--- > Plot Output                                                ----
# Wall
ggplot(outp$wall %>% filter(tokenvalue <= 10), aes(s, LHSG-RHSG, group = airfoil)) +
  geom_line(aes(colour = tokenvalue)) +
  coord_cartesian(ylim = c(-8, 8)) +
  scale_colour_gradientn(colours=rev(spectralpalette(18)), 
                         breaks = seq(2, 20, 2),
                         limits = c(2, 20),
                         name = "N_P") +
  facet_wrap(~airfoilname) +
  ggtitle(expression("Numerical Error on Airfoil Surface, N_P" <= "10"))
ggsave(paste0("NSErr-le10.png"),
       scale = 2, width = 12, height = 6, units = "cm", dpi = 300)

ggplot(outp$wall %>% filter(tokenvalue > 10), aes(s, LHSG-RHSG, group = airfoil)) +
  geom_line(aes(colour = tokenvalue)) +
  coord_cartesian(ylim = c(-8, 8)) +
  scale_colour_gradientn(colours=rev(spectralpalette(18)), 
                         breaks = seq(2, 20, 2),
                         limits = c(2, 20),
                         name = "N_P") +
  facet_wrap(~airfoilname) +
  ggtitle(expression("Numerical Error on Airfoil Surface, N_P" > "10"))
ggsave(paste0("NSErr-ge10.png"),
       scale = 2, width = 12, height = 6, units = "cm", dpi = 300)

ggplot(outp$wall %>% filter(tokenvalue == 10), aes(s, LHSG-RHSG, group = airfoil)) +
  geom_line(aes(colour = tokenvalue)) +
  coord_cartesian(ylim = c(-8, 8)) +
  scale_colour_gradient(low = "black", high = "black", name = "N_P") +
  facet_wrap(~airfoilname) +
  ggtitle("Numerical Error on Airfoil Surface, N_P = 10")
ggsave(paste0("NSErr-eq10.png"),
       scale = 2, width = 12, height = 6, units = "cm", dpi = 300)

# Node
conv <- list()
conv$noder <- outp$node %>% 
  filter(airfoilname == "naca0012r", local != 1) %>% 
  select(u, v, p, dodx, dody, dpdx, dpdy, o, tokenvalue)
conv$nodeo <- outp$node %>% 
  filter(airfoilname == "naca0012", local != 1) %>% 
  select(u, v, p, dodx, dody, dpdx, dpdy, o, tokenvalue)
conv$baser <- filter(conv$noder, tokenvalue == 19) %>% mutate(tokenvalue = 0)
conv$baseo <- filter(conv$nodeo, tokenvalue == 19) %>% mutate(tokenvalue = 0)
conv$noder <- filter(conv$noder, tokenvalue != 19)
conv$nodeo <- filter(conv$nodeo, tokenvalue != 19)

conv$diffr <- (conv$noder - conv$baser[rep(1:nrow(conv$baser), 16),])/
  mutate(conv$baser[rep(1:nrow(conv$baser), 16),], tokenvalue = 1)
conv$diffo <- (conv$nodeo - conv$baseo[rep(1:nrow(conv$baseo), 16),])/
  mutate(conv$baseo[rep(1:nrow(conv$baseo), 16),], tokenvalue = 1)

conv$diffr <- conv$diffr %>% group_by(tokenvalue)
conv$diffo <- conv$diffo %>% group_by(tokenvalue)

conv$diffr[is.na(conv$diffr)] = 0
conv$diffo[is.na(conv$diffo)] = 0

conv$summr <- summarise_all(conv$diffr, funs(max)) %>% mutate(airfoil = "naca0012")
conv$summo <- summarise_all(conv$diffo, funs(max)) %>% mutate(airfoil = "naca0012r")

conv$summary <- rbind(conv$summo, conv$summr) %>% 
  gather(var, norm_inf, -tokenvalue, -airfoil) %>%
  mutate(order = ifelse(var %in% c("u", "v", "p"), "Primary (u, v, p)", "Derivative (o, dpdx, dpdy)")) %>%
  mutate(order = ifelse(var %in% c("dodx", "dody"), "2nd Derivative (dodx, dody)", order))

conv$summary$airfoil <- factor(conv$summary$airfoil, 
                               levels = c("naca0012r", "naca0012"))
conv$summary$var <- factor(conv$summary$var, 
                           levels = c("u", "v", "p", "o", "dpdx", "dpdy", "dodx", "dody"))
conv$summary$order <- factor(conv$summary$order, 
                             levels = c("Primary (u, v, p)", "Derivative (o, dpdx, dpdy)", "2nd Derivative (dodx, dody)"))

ggplot(conv$summary, aes(tokenvalue, norm_inf, group = interaction(var, airfoil), colour = var)) +
  geom_line(aes(linetype = airfoil)) + 
  geom_point(aes(shape = airfoil)) +
  scale_y_continuous(trans='log10') +
  facet_wrap(~order, scales = "free_y") +
  xlab("N_P") +
  scale_linetype_discrete(name = "Airfoil") +
  scale_colour_discrete(name = "Variable") +
  scale_shape_discrete(guide = "none") +
  scale_x_continuous(breaks = seq(4, 20, 4))

ggsave(paste0("node.png"),
       scale = 2.5, width = 12, height = 5, units = "cm", dpi = 300)

ggsave(paste0("node_small.png"),
       scale = 1.8, width = 12, height = 5, units = "cm", dpi = 300)

# BVFa
plot_bvfa <- outp$BVFa %>% 
  select(RHS, LHS, airfoilname, tokenvalue) %>% 
  gather(EqSide, Value, -airfoilname, -tokenvalue)
plot_bvfa$airfoilname <- factor(plot_bvfa$airfoilname, 
                                levels = c("naca0012r", "naca0012"))
plot_bvf <- ggplot(plot_bvfa, aes(tokenvalue, Value,
                      group = interaction(airfoilname, EqSide), colour = airfoilname)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
  geom_path (aes(linetype = airfoilname)) + 
  geom_point(aes(shape = airfoilname)) +
  facet_wrap(~EqSide, scales = "free_y") +
  xlab("N_P Value") + ylab("Integral of BVF")

ggsave(paste0("BVF.png"), plot_bvf,
       scale = 2, width = 12, height = 5, units = "cm", dpi = 300)

# Inte
plot_inte <- outp$inte %>% 
  filter(regn == "Total") %>%
  select(-airfoil, -seshname, -ID, -time, -kinvis, -a) %>% 
  rename(Circulation = o, int.omega.x = ox, int.omega.y = oy) %>% 
  gather(Term, Value, -airfoilname, -tokenvalue, -regn)
plot_inte$airfoilname <- factor(plot_inte$airfoilname, 
                                levels = c("naca0012r", "naca0012"))
ggplot(plot_inte, 
       aes(tokenvalue, Value,
           group = interaction(airfoilname, Term), colour = airfoilname)) +
  geom_path (aes(linetype = airfoilname)) + 
  geom_point(aes(shape = airfoilname)) +
  facet_wrap(~Term, scales = "free_y") +
  xlab("N_P Value") + ylab("Integral of Vorticity")

ggsave(paste0("inte.png"),
       scale = 2, width = 12, height = 5, units = "cm", dpi = 300)


# Junk
ggplot(outp$wall, aes(s, sign(LHSG)*log10(abs(LHSG)), group = airfoil)) +
  geom_line(aes(colour = tokenvalue)) +
  scale_colour_gradientn(colours=rev(spectralpalette(17)), name = "N_P") +
  facet_wrap(~airfoilname)

ggplot(outp$wall, aes(s, LHSG-RHSG, group = airfoil)) +
  geom_line(aes(colour = tokenvalue)) +
  coord_cartesian(ylim = c(-8, 8)) +
  scale_colour_gradientn(colours=rev(spectralpalette(18)), 
                         breaks = seq(2, 20, 2),
                         limits = c(2, 20),
                         name = "N_P") +
  facet_wrap(~airfoilname) +
  ggtitle("Numerical Error on Airfoil Surface")

ggplot(outp$wall %>% filter(tokenvalue <= 10), aes(s, LHSG-RHSG, group = airfoil)) +
  geom_line(aes(colour = tokenvalue)) +
  coord_cartesian(ylim = c(-8, 8)) +
  scale_colour_gradientn(colours=rev(spectralpalette(18)), 
                         breaks = seq(2, 20, 2),
                         limits = c(2, 20),
                         name = "N_P") +
  facet_wrap(~airfoilname) +
  ggtitle(expression("Numerical Error on Airfoil Surface, N_P" <= "10"))
ggplot(outp$wall %>% filter(tokenvalue > 10), aes(s, LHSG-RHSG, group = airfoil)) +
  geom_line(aes(colour = tokenvalue)) +
  coord_cartesian(ylim = c(-8, 8)) +
  scale_colour_gradientn(colours=rev(spectralpalette(18)), 
                         breaks = seq(2, 20, 2),
                         limits = c(2, 20),
                         name = "N_P") +
  facet_wrap(~airfoilname) +
  ggtitle(expression("Numerical Error on Airfoil Surface, N_P" > "10"))
ggplot(outp$wall %>% filter(tokenvalue == 8), aes(s, LHSG-RHSG, group = airfoil)) +
  geom_line(aes(colour = tokenvalue)) +
  coord_cartesian(ylim = c(-8, 8)) +
  scale_colour_gradientn(colours=rev(spectralpalette(18)), 
                         breaks = seq(2, 20, 2),
                         limits = c(2, 20),
                         name = "N_P") +
  facet_wrap(~airfoilname) +
  ggtitle(expression("Numerical Error on Airfoil Surface, N_P" = "8"))
  

ggplot(outp$wall, aes(s, sign(nserrG)*log10(abs(nserrG)), group = airfoil)) +
  geom_line(aes(colour = tokenvalue))


conv <- list()
conv$noder <- outp$node %>% 
  filter(airfoilname == "naca0012r", local != 1) %>% 
  select(u, v, p, dodx, dody, dpdx, dpdy, o, tokenvalue)
conv$nodeo <- outp$node %>% 
  filter(airfoilname == "naca0012", local != 1) %>% 
  select(u, v, p, dodx, dody, dpdx, dpdy, o, tokenvalue)
conv$baser <- filter(conv$noder, tokenvalue == 19) %>% mutate(tokenvalue = 0)
conv$baseo <- filter(conv$nodeo, tokenvalue == 19) %>% mutate(tokenvalue = 0)
conv$noder <- filter(conv$noder, tokenvalue != 19)
conv$nodeo <- filter(conv$nodeo, tokenvalue != 19)

conv$diffr <- (conv$noder - conv$baser[rep(1:nrow(conv$baser), 16),])/
  mutate(conv$baser[rep(1:nrow(conv$baser), 16),], tokenvalue = 1)
conv$diffo <- (conv$nodeo - conv$baseo[rep(1:nrow(conv$baseo), 16),])/
  mutate(conv$baseo[rep(1:nrow(conv$baseo), 16),], tokenvalue = 1)

conv$diffr <- conv$diffr %>% group_by(tokenvalue)
conv$diffo <- conv$diffo %>% group_by(tokenvalue)

conv$diffr[is.na(conv$diffr)] = 0
conv$diffo[is.na(conv$diffo)] = 0

conv$summr <- summarise_all(conv$diffr, funs(max)) %>% mutate(airfoil = "naca0012")
conv$summo <- summarise_all(conv$diffo, funs(max)) %>% mutate(airfoil = "naca0012r")

conv$summary <- rbind(conv$summo, conv$summr) %>% 
  gather(var, norm_inf, -tokenvalue, -airfoil) %>%
  mutate(order = ifelse(var %in% c("u", "v", "p"), "Primary (u, v, p)", "Derivative (o, dpdx, dpdy)")) %>%
  mutate(order = ifelse(var %in% c("dodx", "dody"), "Second Derivative (dodx, dody)", order))

conv$summary$airfoil <- factor(conv$summary$airfoil, 
                               levels = c("naca0012r", "naca0012"))
conv$summary$var <- factor(conv$summary$var, 
                           levels = c("u", "v", "p", "o", "dpdx", "dpdy", "dodx", "dody"))
conv$summary$order <- factor(conv$summary$order, 
                             levels = c("Primary (u, v, p)", "Derivative (o, dpdx, dpdy)", "Second Derivative (dodx, dody)"))

ggplot(conv$summary, aes(tokenvalue, norm_inf, group = interaction(var, airfoil), colour = var)) +
  geom_line(aes(linetype = airfoil)) + 
  geom_point(aes(shape = airfoil)) +
  scale_y_continuous(trans='log10') +
  facet_wrap(~order, scales = "free_y") +
  xlab("N_P") +
  scale_linetype_discrete(name = "Airfoil") +
  scale_colour_discrete(name = "Variable") +
  scale_shape_discrete(guide = "none") +
  scale_x_continuous(breaks = seq(4, 20, 4))
  

