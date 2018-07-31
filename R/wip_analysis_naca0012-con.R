#============================#
# Analysis of NACA0012 Example
# Alwin Wang
#----------------------------#

#--- Set Up                                     ----
# Use rstudioapi to get saved location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
#--- * Scripts                                                    ----
# Source Required Scripts
source("src_library-manager.R")                                 # Call libraries and install missing ones
source("src_helper-functions.R")                                # Smaller functions used
# Additional scripts here
# ggplot2 setup (consider moving to a separate script)
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <-                                              # Custom spectral pallette
  colorRampPalette(rev(brewer.pal(11, "Spectral")))             #  usage: spectralpallette(10) 
#--- * Required Paths                                             ----
if (FALSE) {                                                    # Test cases
  saveplot   = "../src-example/NACA0012_convergence/"             # Pass in later to plot_data
  airfoil    = "NACA0012"
  folderpath = "../src-example/NACA0012_convergence/NACA0012-N_P04/"
  seshpath   = "naca0012-N_P04"
  dumppath   = "naca0012-N_P04-00.dump"
  # Required input dataframes
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
    ID         = "NACA0012-AoA04-N_P5",
    stringsAsFactors = FALSE)
  data_dump <- data.frame(
    data_mesh,
    dumpfile = dumppath,
    dumppath = paste0(folderpath, dumppath),
    stringsAsFactors = FALSE)
  rm(airfoil, folderpath, seshpath, dumppath)
}

#--- Airfoil Calculation                                          ----
# Determine things like spline distance once per unique airfoil i.e. boundary profile and wall output
BatchLoadAirfoil <- function(data_airfoil) {
  source("src_library-manager.R")                                 # Call libraries and install missing ones
  source("src_helper-functions.R")                                # Smaller functions used
  #--- * Boundary Data                                              ----
  # Not sure this is actually used anywhere...
  bndry <- LoadBndry(data_airfoil$folder)
  #--- * Wall Mesh Data                                             ----
  wallmesh <- LoadWallGrad(data_airfoil$seshpath)
  long_wall <- AirfoilLongWall(wallmesh)
  long_wall <- AirfoilSpline(long_wall)
  long_wall <- AirfoilNorm  (long_wall)
  #--- > Airfoil Calc Output                                        ----
  list_airfoil <- list(
    airfoil   = data_airfoil,
    bndry     = bndry,
    long_wall = long_wall)
  rm(data_airfoil, bndry, long_wall, wallmesh)
  # Output
  return(list_airfoil)
}
#--- Session and Mesh Calculation                                 ----
BatchLoadMesh <- function(data_mesh, outp_airfoil) {
  source("src_library-manager.R")                                 # Call libraries and install missing ones
  source("src_helper-functions.R")                                # Smaller functions used
  list_airfoil <- outp_airfoil[[data_mesh$airfoil]]
  #--- * Airfoil Data                                               ----
  long <- list()
  #--- * Session Data                                               ----
  session   <- LoadSeshFileKeywords(data_mesh$seshpath)
  long$sesh <- LongSesh(session)
  #--- * Mesh Data                                                  ----
  long$mesh <- LoadMesh(data_mesh$seshpath)
  long$mesh <- LongMesh(long$mesh, long$sesh)
  #--- * Wall Data                                                  ----
  long$wall <- list_airfoil$long_wall
  long$wall <- LongWall(long$wall, long$mesh)
  #--- * Local Data                                                 ----
  long$mesh <- LocalMesh(long$mesh, long$wall)
  #--- > Sesh & Mesh Calc Output                                    ----
  list_mesh <- list(
    wall = long$wall,
    mesh = long$mesh)
  rm(data_mesh, long, session)
  # output
  return(list_mesh)
}
#--- Dump File Calculation                                        ----
BatchLoadDump <- function(data_dump, outp_mesh) {
  source("src_library-manager.R")                                 # Call libraries and install missing ones
  source("src_helper-functions.R")                                # Smaller functions used
  list_mesh <- outp_mesh[[data_dump$ID]]
  #--- * Dump Data                                                  ----
  dump      <- LoadGradFieldDump(data_dump$folder, data_dump$dumpfile)
  dump$dump <- DumpMesh(list_mesh$mesh, dump$dump)
  dump$wall <- DumpWall(list_mesh$wall, dump$dump)
  #--- * Accleration Data                                           ----
  # Note that tangent direction based on spline calc
  LoadSeshBCEq(data_dump$seshpath, "MOD_ALPHA_X")
  dump$a    <- BC_mod_alpha_x(dump$time)
  dump$wall <- DumpAccel(dump$a, dump$wall)
  #--- * Pressure Data                                              ----
  dump$wall <- DumpPres(dump$wall)
  #--- * Vorticity Data                                             ----
  order = 3 # Ideally, scale should be order/N_P
  dump$offs <- AirfoilOffset(dump$wall, nsteps = order, scale = (order)/data_dump$tokenvalue)
  dump$offs <- DumpVortInterp(dump$offs, dump$dump,
                              linear = TRUE, extrap = FALSE, round = NULL)
  dump$offs <- FiniteDiff(dump$offs, "o", order = order)
  dump$offs <- DumpVortGrad(dump$offs)
  dump$wall <- DumpVortJoin(dump$wall, dump$offs, dump$kinvis)
  # Reduce dump down to just the nodes for field comparison
  dump$node <- filter(dump$dump, node)
  #--- > Dump Calc Output                                           ----
  data_plot <- bind_rows(dump[c("time", "kinvis", "a")])
  data_plot <- cbind(data_dump, data_plot)
  list_dump <- c(
    data_plot = list(data_plot), 
    dump = dump[c("wall", "offs", "node")])
  names(list_dump) <- c("data_plot", "wall", "offs", "node")
  rm(data_dump, data_plot, dump, order)
  # output
  return(list_dump)
}

#--- Batch Calculation                                            ----
# Batch list to process
batchfolder = "../src-example/NACA0012_convergence"
df_batch <- ListSesh(batchfolder)
# List of unique airfoils & N_P
li_airfoil   <- split (df_batch, df_batch$airfoil)
li_airfoil   <- lapply(li_airfoil, function(df) df[1,])
outp_airfoil <- pblapply(li_airfoil, BatchLoadAirfoil)
# List of meshs (should turn this into a function!)
tokenwords = list("N_P")                                        # Find unique bndry and knot N values
df_mesh <- df_batch %>%
  rowwise() %>%
  do(data.frame(., LoadSeshTokenWords(.$seshpath, tokenwords),
                stringsAsFactors = FALSE)) %>%
  mutate(ID = paste0(airfoil))                                  # Determine ID later
li_mesh   <- split (df_mesh, df_mesh$ID)                        # Unique airfoil and N_P
li_mesh   <- lapply(li_mesh, function(x) x[1,])                 # Take only the first row in each
outp_mesh <- pblapply(li_mesh, BatchLoadMesh, outp_airfoil)
# List of dumps
df_dump <- df_mesh %>% 
  rowwise() %>% 
  do(data.frame(., dumpfile = ListDump(.$folder, .$seshname),
                stringsAsFactors = FALSE)) %>%
  mutate(dumppath = paste0(folder, dumpfile))
li_dump   <- split(df_dump, df_dump$dumppath)
outp_dump <- pblapply(li_dump, BatchLoadDump, outp_mesh)
#--- > Batch Calc Output                                        ----
outp <- list()
# Wall
outp$wall <- lapply(outp_dump, function(dump) cbind(dump$wall, 
  select(dump$data_plot, airfoil, seshname, tokenvalue, ID, time, kinvis, a)))
outp$wall <- bind_rows(outp$wall)
# Offs
outp$offs <- lapply(outp_dump, function(dump) cbind(dump$offs, 
  select(dump$data_plot, airfoil, seshname, tokenvalue, ID, time, kinvis, a)))
outp$offs <- bind_rows(outp$offs)
# Node
outp$node <- lapply(outp_dump, function(dump) cbind(dump$node, 
  select(dump$data_plot, airfoil, seshname, tokenvalue, ID, time, kinvis, a)))
outp$node <- bind_rows(outp$node)
rm(df_batch,     df_mesh,   df_dump, 
   li_airfoil,   li_mesh,   li_dump,
   outp_airfoil, outp_mesh, outp_dump,
   tokenwords)

#--- Convergence                                                ----
conv <- list()
# Navier-Stokes for vorticity convergence
conv$wall <- outp$wall %>%
  select(tokenvalue, s, LHSG, RHSG, nserrG) %>%
  group_by(tokenvalue)
conv$nsabserr <- conv$wall %>% 
  summarise(mean = mean(abs(nserrG)),
            max  = max (abs(nserrG)),
            relmean = mean(abs(nserrG/RHSG)),
            relmax  = max (abs(nserrG/RHSG)))

ggplot(conv$wall, 
       aes(s, nserrG, group = tokenvalue, colour = tokenvalue)) +
  geom_point(shape = 'o') +
  geom_line (data = filter(conv$wall, tokenvalue == 10))

ggplot(conv$nsabserr, aes(tokenvalue, mean)) +
  geom_line() + geom_point() +
  scale_y_continuous(trans='log10')
# Convergence of nxG nxyG?
# Convergence of fields
# u v p dodx dody dpdx dpdy o 
conv$fld <- outp$node %>%
  select(tokenvalue, u, v, p, dodx, dody, dpdx, dpdy, o)
maxNP <- max(conv$fld$tokenvalue)
conv_maxN <- filter(conv$fld, tokenvalue == maxNP) %>% 
  mutate(tokenvalue = 0)
conv_err  <- split(conv$fld, conv$fld$tokenvalue)
conv_err  <- lapply(conv_err, function(df) abs(df - conv_maxN))
conv_err  <- bind_rows(conv_err) %>%
  filter(tokenvalue != maxNP) %>% 
  group_by(tokenvalue)
conv$flderr <- summarise_all(conv_err, funs(mean))

conv_flderrplot <- gather(conv$flderr, var, norm_inf, -tokenvalue) %>%
  mutate(order = ifelse(var %in% c("u", "v", "p"), "1", "2")) %>%
  mutate(order = ifelse(var %in% c("dodx", "dody"), "3", order))
ggplot(conv_flderrplot, aes(tokenvalue, norm_inf, group = var, colour = var)) +
  geom_line() + geom_point(aes(shape = order)) +
  scale_y_continuous(trans='log10') +
  facet_wrap(~order, scales = "free_y")
