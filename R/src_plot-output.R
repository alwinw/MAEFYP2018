#============================#
# Plot Output
# Alwin Wang
#----------------------------#

# Custom ggplot2 setup
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

#--- Plot N-S LHS vs RHS ----
PlotNS <- function(dump, dumpval, long, save = TRUE) { # Handle some output variable for path
  # Vertical Lines
  plot_vline <- data.frame(                                       # LE, TE, LE limits for vertical lines
    te_up = min(long$walldata$s),
    le = long$walldata[long$walldata$theta == pi,]$s[1],
    te_lo = max(long$walldata$s))
  plot_surf <- data.frame(                                        # Labels for surfaces
    x = c(plot_vline$te_up,
          (plot_vline$te_up + plot_vline$le)/2,
          plot_vline$le,
          (plot_vline$le + plot_vline$te_lo)/2,
          plot_vline$te_lo),
    y = rep(-20, 5),
    labels = c("TE", "Upper Surface", "LE", "Lower Surface", "TE"))
  plot_xbreaks <-                                                 # Breaks in x axis
    c(seq(plot_vline$te_up, plot_vline$le, length.out = 3)[1:2],
      seq(plot_vline$le, plot_vline$te_lo, length.out = 3))
  # Title
  plot_title <- paste(
    dumpval$airfoil,
    paste("Kinematic Viscosity:", format(dump$kinvis, scientific = FALSE)),
    paste("Time:", sprintf("%4.2f", dump$time)),
    paste("Acceleration:", sprintf("%+6.3f", dump$acceleration)),
    sep = "\n")
  plot_filename <- paste(
    dumpval$ID,
    paste0("v", format(dump$kinvis, scientific = TRUE)),
    paste0("t", sprintf("%05.3f", dump$time)),
    paste0("a", sprintf("%+07.3f", dump$acceleration)),
    "ns",
    sep = "_")
  # Plot
  plot_NS <- ggplot(dump$threaddata %>%                           # Initiate plot with threaddata
                      filter(wall) %>% arrange(s), aes(s)) +
    geom_vline(xintercept = as.numeric(plot_vline),               # Vertical lines for LE, TE, LE
               colour = "grey", linetype = "dashed") +
    geom_label(aes(x, y, label = labels), plot_surf,              # Surface labels
               colour = "grey") +
    geom_path(aes(y = -accels, colour = "dU/dt")) +               # Acceleration terms
    geom_path(aes(y = dpds, colour = "dp/ds")) +                  # Pressure field
    geom_path(aes(y = - accels + dpds, colour = "LHS"),           # LHS, acceleration + pressure
              linetype = "dashed") +
    geom_point(aes(s, t_enum_diff*dump$kinvis, colour = "RHS"),   # RHS, v * dw/dz
               dump$offset, shape = "O") +
    xlab("s") + 
    scale_x_continuous(breaks = plot_xbreaks, 
                       labels = function(x) sprintf("%.2f", x)) +
    ylab(NULL) + ylim(c(-20, 40)) +
    scale_color_manual(
      name = "Legend",
      values = c("dp/ds" = "red", "dU/dt" = "blue", "LHS" = "purple", "RHS" = "purple"),
      labels = c(
        expression(frac(1, rho)~frac(partialdiff*p, partialdiff*s)), 
        expression(frac(partialdiff*U, partialdiff*t)), 
        expression(
          bgroup("(",
                 frac(1, rho)~frac(partialdiff*p, partialdiff*s) +
                   frac(partialdiff*U, partialdiff*t), ")")),
        expression(nu~frac(partialdiff*omega, partialdiff*z))),
      guide = guide_legend(
        override.aes = list(
          linetype = c(rep("solid", 2), "dashed", "blank"),
          shape = c(rep(NA, 3), "O")))) +
    theme(legend.key.size = unit(2.25, "lines"),
          legend.text.align = 0.5,
          legend.direction = "vertical", 
          legend.position = "right",
          legend.background = element_rect(colour = "black", size = 0.3)) +
    ggtitle(plot_title)
  # Save plot
  if (save) {
    plotpath = paste0("../output-plot", "/", dumpval$airfoil, "_", dumpval$seshname)
    if (!dir.exists(plotpath)) dir.create(plotpath, recursive = TRUE)
    ggsave(paste0(plot_filename, ".png"), plot_NS, path = plotpath,
           width = 8, height = 7)
  } else {
    print(plot_NS)
  }
  # Return
  return(NULL)
}

#--- Plot Acceleration ----
PlotAccel <- function(dump, dumpval,save = TRUE) {
  # ASSUME the time is 0 to 2 which may NOT always be the case
  plot_filename <- paste(
    dumpval$ID,
    paste0("v", format(dump$kinvis, scientific = TRUE)),
    paste0("t", sprintf("%05.3f", dump$time)),
    paste0("a", sprintf("%+07.3f", dump$acceleration)),
    "a",
    sep = "_")
  plot_aval <- data.frame(t = seq(0, 2, length.out = 500)) %>%
    mutate(a = BC_mod_alpha_x(t))
  plot_accel <- ggplot(plot_aval, aes(t, a)) +
    geom_path(colour = "blue") +
    geom_point(aes(dump$time, dump$acceleration),
               colour = "blue")
  # Save plot
  if (save) {
    plotpath = paste0("../output-plot", "/", dumpval$airfoil, "_", dumpval$seshname)
    if (!dir.exists(plotpath)) dir.create(plotpath, recursive = TRUE)
    ggsave(paste0(plot_filename, ".png"), plot_accel, path = plotpath,
           width = 6, height = 2)
  } else {
    print(plot_accel)
  }
  # Return
  return(NULL)
}


#--- Plot Leading Edge ----
PlotLE <- function(dump, dumpval, save = TRUE) {
  # Mesh to interpolate over
  n = 250
  plot_le_df = data.frame(
    x = rep(seq(-0.455, -0.245, length.out = n), each = n),
    y = rep(seq(-0.07, 0.1, length.out = n), times = n))
  #--- Relevant enum ----
  # Original mesh
  poly_df <- dump$threaddata %>% 
    filter(x >= -0.8, x <= 0.1,
           y >= -0.3, y <= 0.5,
           seshnode) %>%
    arrange(enum, ncorner) %>%
    select(x, y, enum)
  pts_df <- plot_le_df
  # Split the dataframe into a list based on enum and then remove enum from df in the list
  poly_list <- split(poly_df, poly_df$enum)
  # Convert the list to Polygon, then create a Polygons object
  poly_sp <- sapply(poly_list, function(poly){Polygons(list(Polygon(poly[, c("x", "y")])), ID = poly[1, "enum"])})
  poly_sp <- SpatialPolygons(poly_sp)
  # Convert points to coordinates
  pts_ps <- pts_df
  coordinates(pts_ps) <- ~x+y
  # Determine polygons points are in
  pts_return <- over(pts_ps, poly_sp, returnList = FALSE)
  plot_le_df$enum <- unique(poly_df$enum)[pts_return] 
  plot_le_df <- filter(plot_le_df, !is.na(enum))
  plot_le_list <- split(plot_le_df, plot_le_df$enum)
  #--- Threaddata for enum ----
  plot_full_df <- dump$threaddata %>%
    filter(enum %in% names(plot_le_list)) %>%
    select(x, y, t, enum, nnum, ncorner, local)
  plot_full_list <- split(plot_full_df, plot_full_df$enum)
  
  plot_le_interp <- lapply(names(plot_le_list), function(enum) {
    plot_le_interp <- as.data.frame(
      interpp(x = plot_full_list[[enum]]$x, y = plot_full_list[[enum]]$y, z = plot_full_list[[enum]]$t,
              xo = plot_le_list[[enum]]$x, yo = plot_le_list[[enum]]$y,
              linear = TRUE,
              duplicate = "strip"))
    return(plot_le_interp)
  })
  
  plot_le_interp <- bind_rows(plot_le_interp) %>%
    filter(!is.na(z))
  
  # plot_le_interp <- as.data.frame(
  #   interpp(x = plot_full_df$x, y = plot_full_df$y, z = plot_full_df$t,
  #           xo = plot_le_df$x, yo = plot_le_df$y,
  #           linear = TRUE,
  #           duplicate = "strip"))
  colnames(plot_le_interp) <- c("x", "y", "t")
  
  plot_le <- ggplot(plot_le_interp) +
    geom_tile(aes(x, y, fill = ifelse(abs(t) > 200, 200*sign(t), t))) +
    # stat_contour(aes(x, y, z = t, fill = ..level..), geom = "polygon",
                 # bins = 10) +
    geom_point(aes(x, y, fill = ifelse(abs(t) > 200, 200*sign(t), t), colour = local), 
               data = plot_full_df,
               shape = 21) +
    geom_polygon(aes(x, y, group = enum, colour = local), fill = NA, alpha = 0.1, linetype = "dotted",
                 data = plot_full_df %>% filter(!is.na(nnum)) %>% arrange(enum, ncorner)) +
    geom_polygon(aes(x, y), fill = "white", colour = "black", alpha = 1,
                 data = dump$accel) +
    coord_fixed(
      xlim = c(-0.45, -0.25),
      ylim = c(-0.055, 0.090),
      expand = FALSE) +
    scale_fill_gradientn(name = "vorticity\n", colours = spectralpalette(20), limits = c(-200, 200)) +
    # scale_colour_gradientn(colours = spectralpalette(20), limits = c(-150, 150))
    scale_colour_gradient(low = "grey99", high = "grey80", guide = "none")
  
  plot_filename <- paste(
    dumpval$ID,
    paste0("v", format(dump$kinvis, scientific = TRUE)),
    paste0("t", sprintf("%05.3f", dump$time)),
    paste0("a", sprintf("%+07.3f", dump$acceleration)),
    "le",
    sep = "_")
  
  # Save plot
  if (save) {
    plotpath = paste0("../output-plot", "/", dumpval$airfoil, "_", dumpval$seshname)
    if (!dir.exists(plotpath)) dir.create(plotpath, recursive = TRUE)
    ggsave(paste0(plot_filename, ".png"), plot_le, path = plotpath,
           width = 6, height = 5)
  } else {
    print(plot_le)
  }
  # Return
  return(NULL)
}
