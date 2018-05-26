#============================#
# Plot Output
# Alwin Wang
#----------------------------#

# Custom ggplot2 setup
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

#--- Plot N-S LHS vs RHS ----
PlotNS <- function(dump, dumpval, long, save = TRUE, saveplot = NULL) {
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
    paste("Kinematic Viscosity:", sprintf("%.4f", dump$kinvis)),
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
               dump$offset, shape = "o", alpha = 0.3) +
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
          shape = c(rep(NA, 3), "O"),
          alpha = rep(1, 4)))) +
    theme(legend.key.size = unit(2.25, "lines"),
          legend.text.align = 0.5,
          legend.direction = "vertical", 
          legend.position = "right",
          legend.background = element_rect(colour = "black", size = 0.3)) +
    ggtitle(plot_title)
  # Save plot
  if (save) {
    plotpath = paste0(saveplot, "/", dumpval$airfoil, "_", dumpval$seshname)
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
PlotAccel <- function(dump, dumpval, save = TRUE, saveplot = NULL) {
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
    plotpath = paste0(saveplot, "/", dumpval$airfoil, "_", dumpval$seshname)
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
PlotLE <- function(dump, dumpval, save = TRUE, saveplot = NULL) {
  # Determine correct element for points for interpolation
  n = 250
  # plot_limits <- data.frame(
  #   x = c(-0.55, 0.65),
  #   y = c(-0.3, 0.3))
  plot_limits <- data.frame(
    x = c(-0.455, -0.245),
    y = c(-0.07, 0.1))
  pts_df = data.frame(
    x = rep(seq(plot_limits$x[1], plot_limits$x[2], length.out = n), each = n),
    y = rep(seq(plot_limits$y[1], plot_limits$y[2], length.out = n), times = n))
  poly_df <- dump$threaddata %>% 
    filter(seshnode) %>%
    arrange(enum, ncorner) %>%
    select(x, y, enum)
  pts_df <- PointinElement(pts_df, poly_df)
  
  # Threaddata for enum
  pts_list <- split(pts_df, pts_df$enum)
  plot_full_df <- dump$threaddata %>%
    filter(enum %in% names(pts_list)) %>%
    select(x, y, t, enum, nnum, ncorner, local)
  plot_full_list <- split(plot_full_df, plot_full_df$enum)
  
  plot_le_interp <- lapply(names(pts_list), function(enum) {
    plot_le_interp <- as.data.frame(
      interpp(x = plot_full_list[[enum]]$x, y = plot_full_list[[enum]]$y, z = plot_full_list[[enum]]$t,
              xo = pts_list[[enum]]$x, yo = pts_list[[enum]]$y,
              linear = TRUE,
              duplicate = "strip"))
    return(plot_le_interp)
  })
  
  plot_le_interp <- bind_rows(plot_le_interp) %>%
    filter(!is.na(z))
  colnames(plot_le_interp) <- c("x", "y", "t")
  
  PlotScaleLimit <- function(t, max) ifelse(abs(t) > max, max*sign(t), t)
  PlotContourBreak <- function(max, step) seq(-max, max, by = step)[seq(-max, max, by = step) != 0]
  
  plot_le <- ggplot(plot_le_interp) +
    # Vorticity tile
    geom_tile(aes(x, y, fill = PlotScaleLimit(t, 200)), alpha = 1) +
    # Original Points
    # geom_point(aes(x, y, fill = PlotScaleLimit(t, 200)), 
    # shape = 23, colour = "white",
    # data = plot_full_df) +
    # Original Mesh
    geom_polygon(aes(x, y, group = enum), 
                 fill = NA, colour = "grey90", alpha = 0.1, linetype = "1616",
                 data = plot_full_df %>% filter(!is.na(nnum)) %>% arrange(enum, ncorner)) +
    # Vorticity contours
    # stat_contour(aes(x, y, z = t, fill = ..level..), geom = "polygon", alpha =1,
    # breaks = PlotContourBreak(200, 10)) +
    stat_contour(aes(x, y, z = t, colour = ..level..), alpha = 1,
                 breaks = PlotContourBreak(200, 10)) +
    # stat_contour(aes(x, y, z = t, colour = -abs(..level..)),
                 # breaks = PlotContourBreak(200, 10)) +
    # Chosen offset points
    # geom_point(aes(x, y, fill = PlotScaleLimit(t_enum, 200)),
    #            data = dump$offset,
    #            shape = 21) +
    # Airfoil
    geom_polygon(aes(x, y), fill = "white", colour = "black", alpha = 1,
                 data = dump$accel) +
    coord_fixed(
      xlim = plot_limits$x,
      ylim = plot_limits$y,
      expand = FALSE) +
    scale_fill_gradientn(name = "vorticity\n", colours = spectralpalette(20), limits = c(-200, 200)) +
    scale_colour_gradientn(colours = spectralpalette(20), limits = c(-200, 200), guide = "none")
    # scale_colour_continuous(low = "grey90", high = "white", guide = "none")
  
  plot_filename <- paste(
    dumpval$ID,
    paste0("v", format(dump$kinvis, scientific = TRUE)),
    paste0("t", sprintf("%05.3f", dump$time)),
    paste0("a", sprintf("%+07.3f", dump$acceleration)),
    "le",
    sep = "_")
  
  # Save plot
  if (save) {
    plotpath = paste0(saveplot, "/", dumpval$airfoil, "_", dumpval$seshname)
    if (!dir.exists(plotpath)) dir.create(plotpath, recursive = TRUE)
    ggsave(paste0(plot_filename, ".png"), plot_le, path = plotpath,
           width = 6, height = 5)
  } else {
    print(plot_le)
  }
  # Return
  return(NULL)
}


#--- Plot LHS RHS Error ----
PlotError <- function(dump, surface, dumpval, long, save = TRUE, saveplot = NULL) {
  summary <- data.frame(
    x = c(2.5, 3),
    y = rep(-90, 2),
    labels = c(paste("median:", sprintf("%.2f%%", median(filter(surface, up)$error)), "\n",
                     "iqr:", sprintf("%.2f%%", IQR(filter(surface, up)$error))),
               paste("median:", sprintf("%.2f%%", median(filter(surface, !up)$error)), "\n",
                     "iqr:", sprintf("%.2f%%", IQR(filter(surface, !up)$error)))))
  
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
    y = rep(-100, 5),
    labels = c("TE", "Upper Surface", "LE", "Lower Surface", "TE"))
  plot_xbreaks <-                                                 # Breaks in x axis
    c(seq(plot_vline$te_up, plot_vline$le, length.out = 3)[1:2],
      seq(plot_vline$le, plot_vline$te_lo, length.out = 3))
  plot_xlabs <- sprintf("%.2f", plot_xbreaks)
  plot_xbreaks <- c(plot_xbreaks, 2.5, 3)
  plot_xlabs <- c(plot_xlabs, "Uppper", "Lower")
  
  ploterror <- ggplot(surface) +
    geom_vline(xintercept = as.numeric(plot_vline),               # Vertical lines for LE, TE, LE
               colour = "grey", linetype = "dashed") +
    geom_col(aes(s, error, fill = up)) +
    geom_boxplot(aes(as.numeric(up)*-0.5 + 3, error, colour = up),
                 outlier.shape = "x") +
    geom_label(aes(x, y, label = labels), plot_surf,              # Surface labels
               colour = "grey40", size = 3) +
    geom_label(aes(x, y, label = labels), summary,                # Summary labels
               colour = "grey40", size = 3) +
    scale_x_continuous(breaks = plot_xbreaks, 
                       labels = plot_xlabs) +
    coord_cartesian(ylim = c(-100, 100)) +
    ylab(expression(frac(RHS - LHS, RHS)%*%100~"%")) +
    guides(fill = FALSE, colour = FALSE) +
    theme(
      axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 0, r = 10000, b = 0, l = 0))
      )
  
  plot_filename <- paste(
    dumpval$ID,
    paste0("v", format(dump$kinvis, scientific = TRUE)),
    paste0("t", sprintf("%05.3f", dump$time)),
    paste0("a", sprintf("%+07.3f", dump$acceleration)),
    "err",
    sep = "_")
  
  # Save plot
  if (save) {
    plotpath = paste0(saveplot, "/", dumpval$airfoil, "_", dumpval$seshname)
    if (!dir.exists(plotpath)) dir.create(plotpath, recursive = TRUE)
    ggsave(paste0(plot_filename, ".png"), ploterror, path = plotpath,
           width = 8, height = 3)
  } else {
    print(ploterror)
  }
  # Return
  return(NULL)
}
