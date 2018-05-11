#============================#
# Vorticity Generation
# Alwin Wang
#----------------------------#

#--- Convert session into a session long mesh ----
LongSession <-function(session) {
  # Elements
  # Check that all elements are quadrangles
  if (as.logical(sum(session$elements$shapetag != "<Q>"))) warning("Not all elements are Quadrangles") 
  # Manipulate data
  long_seshdata <- session$elements %>%
    # Remove shape tag since checked shape already
    select(-shapetag) %>%
    # Gather nodes 1-5 into a single column
    gather(ncorner, nnum, -enum) %>%
    # Combine with nodes to get (x,y) for each element corner
    left_join(., session$nodes, "nnum") %>%
    # Order data
    group_by(enum) %>%
    # Add columns for labels
    mutate(
      elabx = mean(x),
      elaby = mean(y)) %>%
    # Close the loop (without effecting mean)
    rbind(., mutate(filter(., ncorner=="n1"), ncorner = "n5")) %>%
    arrange(enum, ncorner) %>% 
    # Calculate area
    mutate(area = x*lead(y) - lead(x)*y) %>%
    mutate(area = 1/2*abs(sum(ifelse(is.na(area), 0, area)))) %>%
    # Remove extra ncorner
    filter(ncorner != "n5") %>%
    ungroup()
  # Return the output
  return(long_seshdata)
}

#--- Convert mesh into a long mesh ----
LongMesh <- function(mesh, airfoildata) {
  # Determine which mesh data belong to which
  mesh$mnum = 1:nrow(mesh)
  long_meshdata = mesh
  # Return the output
  return(long_meshdata)
}

#--- Determine local mesh ----
LocalWallH <- function(long) {
  # Average height of offset
  long_ave <- long$threaddata %>% 
    group_by(enum) %>%                                              # Group by sesh element
    filter(!is.na(ncorner)) %>%                                     # Only corners of sesh elements
    add_count(wall) %>%                                             # Count number of corners on wall
    filter(n == 2, wall) %>%                                        # Only elements with edge (2 corners) on wall
    mutate(base = EucDist(x, y),                                    # Base length
           aveh = area/base,                                        # Average height
           ar = base/aveh) %>%                                      # Aspect ratio of element
    filter(!is.na(base)) %>%                                        # Reduce to one row per element
    ungroup()
  # Summarise
  summary <- data.frame(
    min = min(long_ave$aveh), 
    median = median(long_ave$aveh),
    mean = mean(long_ave$aveh))
  # Return
  return(c(long, list(wallaveh = long_ave, wallavesum = summary)))
}

#--- Join airfoil data with main data ----
LongAirfoil <- function(long) {
  # Join in initial wall and up varaibles
  long$threaddata = LongJoin(
    long$threaddata, 
    select(long$walldata, x, y, wall, up))
  
  check <- long$threaddata %>% 
    filter(!is.na(ncorner)) %>%
    group_by(enum) %>%
    mutate(
      wall = ifelse(is.na(wall), FALSE, wall),
      checkwall = ifelse(is.na(wall), 0, sum(wall)),
      checkup = sum(up)) 
  
  check_nnum <- unique(check[check$checkwall == 1 & check$wall,]$nnum)
  if (length(check_nnum) != 1) warning("Multiple TE points found")
  check_enum <- check[check$nnum %in% check_nnum,]$enum
  
  checkte <- check[check$enum %in% check_enum,] %>%
    filter(!is.na(ncorner)) %>% 
    group_by(enum) %>% arrange(enum) %>%
    mutate(numup = sum(up, na.rm = TRUE)) %>%
    # Fix up lower element on surface (336)
    mutate(up = ifelse(up == TRUE & checkwall == 2 & numup == 1, FALSE, up)) %>%
    mutate(numup = sum(up, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(nnum) %>% arrange(nnum) %>%
    mutate(eleup = ifelse(nnum == check_nnum, NA, mean(numup))) %>%
    ungroup() %>%
    group_by(enum) %>% arrange(enum) %>%
    mutate(eleup = sum(eleup - 1, na.rm = TRUE)) %>%
    # Fix up lower element behind trailing edge
    mutate(up = ifelse(up == TRUE & eleup == -0.5, FALSE, up)) %>%
    # Clean up
    select(enum, ncorner, nnum, up) %>%
    rename(fixed_up = up)
  
  long_fixte <- LongJoin(long$threaddata, checkte) %>%
    mutate(up = ifelse(!is.na(fixed_up), fixed_up, up)) %>%
    select(-fixed_up)
}


LocalMesh <- function(long) {
  # First nodes and elements
  localmesh <- filter(long$threaddata, wall, !is.na(ncorner))
  localmesh <- list(
    nodes = unique(localmesh$nnum),
    elements = unique(localmesh$enum))
  localmesh$mesh = long$seshdata %>% select(enum, nnum)
  localmesh$mesh$local = 0
  localmesh$mesh$local = 
    ifelse(localmesh$mesh$enum %in% localmesh$elements, 
           localmesh$mesh$local - 1, 
           localmesh$mesh$local)
  # Loop
  while (sum(localmesh$mesh$local == 0) > 0) {
    localmesh$nodes = unique(
      localmesh$mesh[localmesh$mesh$enum %in% localmesh$elements,]$nnum)
    localmesh$elements = unique(
      localmesh$mesh[localmesh$mesh$nnum %in% localmesh$nodes,]$enum)
    localmesh$mesh$local = 
      ifelse(localmesh$mesh$enum %in% localmesh$elements, 
             localmesh$mesh$local - 1, 
             localmesh$mesh$local)
  }
  localmesh$mesh$local = localmesh$mesh$local - min(localmesh$mesh$local) + 1
  # Return
  return(localmesh$mesh)
}
