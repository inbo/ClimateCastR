#' Data thining function
#'
#' A function to reduce the number of data points in a data frame by proximity.
#' The function groups points that are close to each other. The group is then
#' represented by the mean of the points in the group. Closeness is defined by a
#' spatial cutoff distance.
#'
#' @param data A sf object resulting from previous steps in the ClimateCastR - flow
#' @param cutoff A numeric value defining the distance cutoff for grouping points
#' @param fun A function to apply to the grouped points. Default is 'mean'
#'
#' @details
#' The cutoff is set to 1000 meters by default. This means that points that are
#' within 1000 meters of each other will be grouped together. The cutoff can be
#' changed by setting the cutoff argument to a different value.
#'
#' The fun argument can be any function that takes a numeric vector as input and
#' returns a single numeric value. The default is the mean function. Other options
#' are median, min, max & sum.
#'
#' @return A sf object with the same columns as the input data, but with fewer rows
#'
#' @examples
#' \dontrun{
#' # get data using get_data functions or as a result of the data_prep function.
#' taxon_key <- c(2865504, 5274858)
#' df<-get_gbif_data(taxon_key)
#' }
#'

data_thin <- function(df,
                      cutoff = 1000,
                      fun = 'mean') {
  # Check if the input data is an sf object
  if (!inherits(df, "sf")) {
    stop("The input data must be an sf object")
  }

  # Check if the input data has a geometry column
  if (!"geometry" %in% names(df)) {
    stop("The input data must have a geometry column")
  }

  # Check if the input data has a taxonkey column
  if (!"acceptedTaxonKey" %in% names(df)) {
    stop("The input data must have a acceptedTaxonKey column")
  }

  # Check if the input data has a n_obs column
  if (!"n_obs" %in% names(df)) {
    stop("The input data must have a n_obs column")
  }

  # Check if the input data has a decimalLongitude column
  if (!"decimalLongitude" %in% names(df)) {
    stop("The input data must have a decimalLongitude column")
  }

  # Check if the input data has a decimalLatitude column
  if (!"decimalLatitude" %in% names(df)) {
    stop("The input data must have a decimalLatitude column")
  }

  # Check if fun is a accepted function
  if (!fun %in% c("mean", "median", "min", "max", "sum")) {
    stop("The fun argument must be one of 'mean', 'median', 'min', 'max', or 'sum'")
  }

  # get unique taxonkeys
  taxonkeys <- unique(df$acceptedTaxonKey)

  # for loop to go through each taxonkey
  for (t in 1:length(taxonkeys)){
    # subset the data for each taxonkey
    taxa_data <- df %>%
      dplyr::filter(acceptedTaxonKey == taxonkeys[t])

    # get unique year_cat
    year_cats <- unique(taxa_data$year_cat)

    # for loop to go through each year_cat
    for (y in 1:length(year_cats)){
      # subset the data for each year_cat
      data <- taxa_data %>%
        dplyr::filter(year_cat == year_cats[y])

      # Check if the input data has more than one point
      if (nrow(data) < 2) {
        warning("The input data has fewer than two points for taxonkey ", taxonkeys[t], " and year_cat ", year_cats[y])
        data_thin <- data
        next
      }

      # Calculate pairwise distances between points
      dists <- sf::st_distance(data) %>%
        units::drop_units()

      # Set diagonal to NA
      diag(dists) <- NA

      # Create a logical matrix where distances are less than or equal to cutoff
      close <- dists <= cutoff

      # Create a vector to store the group number for each point
      group <- rep(NA, nrow(data))

      # Initialize the group number
      g <- 1

      # Loop through the rows of the logical matrix
      for (i in 1:nrow(close)) {
        # If the point is not already in a group
        if (is.na(group[i])) {
          # Assign the group number to the point
          group[i] <- g
          # Find the points that are close to the current point
          close_points <- which(close[i, ])
          # Assign the group number to the close points
          group[close_points] <- g
          # Increment the group number
          g <- g + 1
        }
      }

      # Create a data frame with the group numbers
      groups <- data.frame(group = group)

      # Calculate the mean of the points in each group
      means <- data %>%
        sf::st_drop_geometry() %>%
        dplyr::mutate(group = groups$group) %>%
        dplyr::group_by(group) %>%
        dplyr::summarize(mean_decimalLongitude = mean(decimalLongitude, na.rm = TRUE),
                         mean_decimalLatitude = mean(decimalLatitude, na.rm = TRUE)) %>%
        sf::st_as_sf(coords = c("mean_decimalLongitude", "mean_decimalLatitude"),
                     crs = 4326,
                     remove = FALSE)

      # add groups to data
      data$group <- groups$group

      # Merge the means with the original data
      data_thin <- dplyr::left_join(data %>% sf::st_drop_geometry(),
                                    means,
                                    by = "group")
      # rbind data for each year_cat
      if(y == 1){
        data_thin_all_year <- data_thin
      } else {
        data_thin_all_year <- rbind(data_thin_all_year, data_thin)
      }
    }
    # rbind data for each taxonkey
    if(t == 1){
      data_thin_all <- data_thin_all_year
    } else {
      data_thin_all <- rbind(data_thin_all, data_thin_all_year)
    }
  }

  if(fun == "mean"){
    dplyr::n_distinct(data_thin_all$group, data_thin_all$acceptedTaxonKey, data_thin_all$year_cat)

    data_thin_redux <- data_thin_all %>%
      dplyr::group_by(group, acceptedTaxonKey, year_cat) %>%
      dplyr::summarise(n_obs = mean(n_obs, na.rm = TRUE),
                       coordinateUncertaintyInMeters = mean(coordinateUncertaintyInMeters, na.rm = TRUE),
                       decimalLatitude = mean_decimalLatitude,
                       decimalLongitude = mean_decimalLongitude,
                       .groups = "keep") %>%
      dplyr::ungroup()
  }
  return(data_thin_all)
}


