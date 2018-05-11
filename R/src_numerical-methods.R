#============================#
# Numerical Methods
# Alwin Wang
#----------------------------#

#--- Required Functions ----
# Heavyside function (step function)
heav <- function(t) ifelse(t>0,1,0)

# Distance function
EucDist <- function(x, y) sqrt((x - lag(x))^2 + (y - lag(y))^2)

# Unique values
UniLeftJoin <- function(longdata, unicols = c("x", "y")) {
  longdata[!duplicated(longdata[,unicols]),]
}

# Custom left join
LongJoin <- function(left_data, right_data, unicols = NULL, wall = FALSE) {
  # Prepare data for joining
  join_cols <-                                                    # Determine columns to join over
    colnames(left_data)[colnames(left_data) %in% colnames(right_data)]
  if (is.null(unicols)) unicols = join_cols                       # Cols to determine duplicate rows
  uni_right_data <- UniLeftJoin(right_data, unicols)              # Remove duplicate rows
  # Join left_data and right_data
  join_data <- left_join(                                         # Join ALL columns and rows
    left_data, uni_right_data, by = join_cols) 
  if (wall) join_data <- mutate(join_data, wall = !is.na(wnum))   # Update wall column if wanted
  # Check for duplicates
  if(nrow(join_data) != nrow(left_data)) {                        # Check the number of rows has not increased
    warning(paste(
      "join_data", "has more rows than", 
      deparse(substitute(left_data)),
      "by", nrow(join_data) - nrow(left_data)))}
  # Check for missing values
  join_anti <- anti_join(                                         # Find rows not merged in
    uni_right_data, left_data, by = join_cols)
  if(nrow(join_anti) != 0) {
    warning(paste(
      deparse(substitute(right_data)), 
      "missing rows in", "join_data",
      "by", nrow(join_anti)))}
  # Return Result
  return(join_data)
}

#--- Cubic Spline Calculus ----
# Determine derivatives and antiderivatives of cubic splines
CubicSplineCalc <- function(cs, order = 0) {
  # Object to return
  ccs <- cs
  # nth Derivative or nth Integral
  while(order != 0) {
    if(order < 0) {
      # Derivative
      if (ncol(ccs$coefs) > 2) ccs$coefs = t(apply(ccs$coefs, 1, polyder))
      else if (ncol(ccs$coefs) == 2) ccs$coefs = apply(ccs$coefs, 1, polyder)
      else ccs$coefs = rep(0, length(ccs$coefs))
      # Ensure correct data structure
      ccs$coefs <- as.matrix(ccs$coefs)
      # Increase order
      order = order + 1}
    else {
      # Integrate
      ccs$coefs = t(apply(ccs$coefs, 1, polyint))
      # Set correct intercepts
      int = c(0)
      for (i in 2:(nrow(ccs$coefs))) {
        int[i] = polyval(ccs$coefs[i-1,], ccs$breaks[i] - ccs$breaks[i-1])
      }
      int = cumsum(int)
      ccs$coefs[,ncol(ccs$coefs)] = int
      # Decrease order
      order = order -1}
  }
  # Adjust to correct order
  ccs$order = ncol(ccs$coefs)
  # Plots
  # plot <- data.frame(x = seq(cs$breaks[1], cs$breaks[length(cs$breaks)], length = 100))
  # plot$y = ppval(cs, plot$x)
  # plot$cy = ppval(ccs, plot$x)
  # ggplot(plot, aes(x = x)) +
  #   geom_path(aes(y = y)) + geom_point(aes(y = y)) +
  #   geom_path(aes(y = cy)) + geom_point(aes(y = cy))
  return(ccs)
}

#--- Interpolation Function ----
# omesh = long_localdump; imesh = airfoildata$offset; onames = c("x", "y", "t", "wall"); inames = c("x", "y"); wallsplit = TRUE
interpolate <- function(omesh, imesh,
                        onames = c("x", "y", "z", "wall"), inames = c("x", "y", "n"),
                        wallsplit = TRUE) {
  # Selet data based on onames for omesh
  omesh <- omesh[onames] 
  colnames(omesh) <- c("x", "y", "z", "wall")
  # rename columns for inames 
  
  # Check for NAs <-- what did I mean for this initially? check imesh = omesh on wall?
  # Split datat into wall and non-wall
  if (wallsplit) {
    # Split the interp mesh into two separate meshes
    wmesh <- filter(imesh, wall)
    imesh <- filter(imesh, !wall)
  }
  # Interpolate onto imesh using duplicate = "strip"
  imesho <- as.data.frame(
    interpp(x = omesh$x, y = omesh$y, z = omesh$z,
            xo = imesh$x, yo = imesh$y,
            linear = FALSE,
            duplicate = "strip")) %>%
    # Make more robust later
    cbind(., imesh[-(1:2)]) %>%
    arrange(snum, nstep)
  # Check for NA
  if (sum(is.na(imesho$z)) > 0) warning("NA found in interpolation")
  # Recombine wall and non wall if wallsplit == TRUE
  if (wallsplit) {
    # Determine wall values for wmesh
    # semi_join(omesh, wmesh, by = c("x", "y"))
    wmesh <- left_join(
      wmesh, omesh[!duplicated(omesh[c("x", "y")]),], 
      by = c("x", "y", "wall"))
    # Join with previous imesho
    imesho <- rbind(imesho, wmesh) %>%
      arrange(snum, nstep)
  }
  # Plot
  ggplot(imesho, aes(x, y, colour = z)) +
    geom_point() + coord_fixed()
  # Check the interpolation accuracy by interpolating back onto omesh and save column % error
  omeshc <- as.data.frame(
    interpp(x = imesho$x, y = imesho$y, z = imesho$z,
            xo = omesh$x, yo = omesh$y,
            linear = FALSE,
            duplicate = "strip")) %>%
    rename(zo = z) %>%
    cbind(., omesh[c("z", "wall")]) %>%
    filter(!is.na(zo)) %>%
    mutate(zerror = (zo - z),
           zrelerror = zerror/z)
  ggplot(omeshc) + geom_density(aes(zerror)) + facet_wrap(~wall, scales = "free_y")
  ggplot(omeshc, aes(x, zerror)) + geom_point() + facet_wrap(~wall, scales = "free_y")
  ggplot(omeshc, aes(x, y, colour = zerror)) +
    geom_point() + coord_fixed()
  # Add original column names
  
  # Check no increase in row number for omesh
  return(imesho)
}