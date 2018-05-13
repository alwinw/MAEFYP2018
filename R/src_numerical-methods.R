#============================#
# Numerical Methods
# Alwin Wang
#----------------------------#

#--- Required Functions ----
# Heavyside function (step function)
heav <- function(t) ifelse(t>0,1,0)

# Distance function
EucDist <- function(x, y) sqrt((x - lag(x))^2 + (y - lag(y))^2)

# Function to find minimum distance
MinS <- function(x, y, csx, csy, dcsx, dcsy, lim_low, lim_up) {
  # Determine minimum distance (bounded!)
  min.out <- fminbnd(
    function(s) {sqrt((x - ppval(csx, s))^2 + (y - ppval(csy, s))^2)},
    lim_low, lim_up)
  # Determine vectors
  vecpt = data.frame(
    x = x - ppval(csx, min.out$xmin),
    y = y - ppval(csy, min.out$xmin))
  vecsf = data.frame(
    x = ppval(dcsx, min.out$xmin),
    y = ppval(dcsy, min.out$xmin))
  # Determin dot products
  dotprod = vecpt$x*vecsf$x + vecpt$y*vecsf$y
  crossprod = vecpt$x*vecsf$y - vecsf$x*vecpt$y
  dist = sqrt(vecpt$x^2 + vecpt$y^2) * sqrt(vecsf$x^2 + vecsf$y^2)
  # Note: if I do dotprod/dist, at the surface ~ 0/0 != 0 because of machine error
  # Create and return output
  return(data.frame(
    stream = min.out$xmin,
    norm = min.out$fmin,
    prec = min.out$estim.prec,
    dotprod = dotprod,
    crossprod = crossprod,
    dist = dist))
}

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

#--- Finite Difference Method ----
FiniteDiff <- function(offset, var) {
  
}



