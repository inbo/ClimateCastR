#' Data thining function
#'
#' A function to reduce the number of data points in a data frame by proximity.
#' The function groups points that are close to each other. The group is then
#' represented by the mean of the points in the group. Closeness is defined by a
#' spatial cutoff distance.
#'
#' @param x A sf object resulting from previous steps in the ClimateCastR - flow
#' @param cutoff A numeric value defining the distance cutoff for grouping points
#' @param fun A function to apply to the grouped points. Default is 'sum'
#' @param n An integer defining the number of points to process at a time. Default is 2000
#'
#' @details
#' The cutoff is set to 1000 meters by default. This means that points that are
#' within 1000 meters of each other will be grouped together. The cutoff can be
#' changed by setting the cutoff argument to a different value.
#'
#' The fun argument can be any function that takes a numeric vector as input and
#' returns a single numeric value. The default is the sum function. Other options
#' are median, min, max, n & mean.
#'
#' The n argument is used to split the data into chunks. This is useful when the
#' data *defined by a taxonkey and year_cat* has more than n points. The default is 2000.
#' When the data has less than n points, the function will process all the data at once.
#' Spatial groups are calculated within a chunk. The function will then merge the
#' chunks and calculate the mean of the points in each group.
#' Increasing the n argument will slow down the function, but may thin the data more.
#' Decreasing the n argument will speed up the function, but may thin the data less.
#'
#' @return A sf object with the same columns as the input data, but with fewer rows
#'
#' @examples
#' \dontrun{
#' # get data using get_data functions or as a result of the data_prep function.
#' taxon_key <- c(2865504, 5274858)
#' df<-get_gbif_data(taxon_key)
#'
#' # thin the data using default values
#' df_thinned <- data_thin(df)
#'
#' # thin the data using a different cutoff value
#' df_thinned <- data_thin(df, cutoff = 500)
#'
#' # thin the data using a different function
#' df_thinned <- data_thin(df, fun = 'mean')
#'
#' # thin the data using a different n value
#' df_thinned <- data_thin(df, n = 1000)
#' }
#'

data_thin <- function(x,
                      cutoff = 1000,
                      fun = 'sum',
                      n = 2000) {
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
  if (!fun %in% c("mean", "median", "min", "max", "sum", "n")) {
    stop("The fun argument must be one of 'mean', 'median', 'min', 'max', 'n' or 'sum'")
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
      message(paste0("Processing taxonkey ", taxonkeys[t], " in period ", year_cats[y]))
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
      message("Calculating pairwise distances between points")
      if(nrow(data) > n){
        message("The data has more than ", n, " points. This may take a while.")
        # split the data into chunks

        chunks <- split(data, 1:nrow(data) %/% n)

        groups <- pbapply::pblapply(seq_along(chunks), function(i) {
          # determine the chunk_number & subset the data
          x <- chunks[[i]]
          chunk_number <- i

          # Calculate pairwise distances between points
          dist_matrix <- sf::st_distance(x) %>%
            units::drop_units()

          # Set diagonal to NA
          diag(dist_matrix) <- NA

          # Create a logical matrix where distances are less than or equal to cutoff
          close <- dist_matrix <= cutoff

          # Create a vector to store the group number for each point
          group <- rep(NA, nrow(x))

          # Initialize the group number
          g <- 1

          # Loop through the rows of the logical matrix
          for (i in 1:nrow(close)) {
            # If the point is not already in a group
            if (is.na(group[i])) {
              # Assign the group number to the point
              group[i] <- paste0(chunk_number, "_", g)
              # Find the points that are close to the current point
              close_points <- which(close[i, ])
              # Assign the group number to the close points
              group[close_points] <- paste0(chunk_number, "_", g)
              # Increment the group number
              g <- g + 1
            }
          }
          return(data.frame(group = group))
        }) %>%
          do.call(rbind, .)

      } else {

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
      }

      # Calculate the mean of the points in each group
      means <- data %>%
        sf::st_drop_geometry() %>%
        dplyr::mutate(group = groups$group) %>%
        dplyr::group_by(group) %>%
        dplyr::summarize(mean_decimalLongitude = mean(decimalLongitude, na.rm = TRUE),
                         mean_decimalLatitude = mean(decimalLatitude, na.rm = TRUE)) %>%
        sf::st_as_sf(coords = c("mean_decimalLongitude", "mean_decimalLatitude"),
                     crs = 4326,
                     remove = FALSE) %>%
        dplyr::ungroup()

      # add groups to data
      data$group <- groups$group

      # Merge the means with the original data
      data_thin <- dplyr::left_join(data %>% sf::st_drop_geometry(),
                                    means,
                                    by = "group")

      n_groups <- dplyr::n_distinct(data_thin$group, data_thin$acceptedTaxonKey, data_thin$year_cat)

      if(n_groups == 1 & nrow(data_thin) > 1){
        warning(paste0("The data for ", taxonkeys[t], " in period ", year_cats[y], " has only one group. Try decreasing the cutoff value."))
      }
      if(n_groups == 0){
        stop("The data has no groups. Data thining was unsuccessful. Try decreasing the cutoff value.")
      }

      if(n_groups < nrow(data_thin)){
        message("Merging groups")
        groups <- unique(data_thin$group)
        # for loop to go through each group
        # initiate progress bar
        pb <- progress::progress_bar$new(format = "  [:bar] :percent ETA: :eta",
                                         total = n_groups,
                                         clear = FALSE,
                                         width = 60)

        for (g in groups){
          pb$tick()
          data_thin_group <- data_thin %>%
            dplyr::filter(group == g)

          if(nrow(data_thin_group) > 1){
            if(fun == "mean"){
              mean_n_obs <- mean(data_thin_group$n_obs, na.rm = TRUE)

              year_cat <- year_cats[y]
              taxonkey <- taxonkeys[t]

              data_thin <- data_thin %>%
                dplyr::filter(group != g) %>%
                dplyr::add_row(
                  year_cat = year_cat,
                  acceptedTaxonKey = taxonkey,
                  n_obs = mean_n_obs,
                  coordinateUncertaintyInMeters = cutoff,
                  decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                  decimalLongitude = data_thin_group$mean_decimalLongitude[1],
                  group = g,
                  geometry = data_thin_group$geometry[1],
                  acceptedScientificName = data_thin_group$acceptedScientificName[1],
                  mean_decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                  mean_decimalLongitude = data_thin_group$mean_decimalLongitude[1])
            }
            if(fun == "median"){
              median_n_obs <- median(data_thin_group$n_obs, na.rm = TRUE)

              year_cat <- year_cats[y]
              taxonkey <- taxonkeys[t]

              data_thin <- data_thin %>%
                dplyr::filter(group != g) %>%
                dplyr::add_row(year_cat = year_cat,
                               acceptedTaxonKey = taxonkey,
                               n_obs = median_n_obs,
                               coordinateUncertaintyInMeters = cutoff,
                               decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                               decimalLongitude = data_thin_group$mean_decimalLongitude[1],
                               group = g,
                               geometry = data_thin_group$geometry[1],
                               acceptedScientificName = data_thin_group$acceptedScientificName[1],
                               mean_decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                               mean_decimalLongitude = data_thin_group$mean_decimalLongitude[1])
            }
            if(fun == "min"){
              min_n_obs <- min(data_thin_group$n_obs, na.rm = TRUE)

              year_cat <- year_cats[y]
              taxonkey <- taxonkeys[t]

              data_thin <- data_thin %>%
                dplyr::filter(group != g) %>%
                dplyr::add_row(year_cat = year_cat,
                               acceptedTaxonKey = taxonkey,
                               n_obs = min_n_obs,
                               coordinateUncertaintyInMeters = cutoff,
                               decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                               decimalLongitude = data_thin_group$mean_decimalLongitude[1],
                               group = g,
                               geometry = data_thin_group$geometry[1],
                               acceptedScientificName = data_thin_group$acceptedScientificName[1],
                               mean_decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                               mean_decimalLongitude = data_thin_group$mean_decimalLongitude[1])
            }
            if(fun == "max"){
              max_n_obs <- max(data_thin_group$n_obs, na.rm = TRUE)

              year_cat <- year_cats[y]
              taxonkey <- taxonkeys[t]

              data_thin <- data_thin %>%
                dplyr::filter(group != g) %>%
                dplyr::add_row(year_cat = year_cat,
                               acceptedTaxonKey = taxonkey,
                               n_obs = max_n_obs,
                               coordinateUncertaintyInMeters = cutoff,
                               decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                               decimalLongitude = data_thin_group$mean_decimalLongitude[1],
                               group = g,
                               geometry = data_thin_group$geometry[1],
                               acceptedScientificName = data_thin_group$acceptedScientificName[1],
                               mean_decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                               mean_decimalLongitude = data_thin_group$mean_decimalLongitude[1])
            }
            if(fun == "sum"){
              sum_n_obs <- sum(data_thin_group$n_obs, na.rm = TRUE)

              year_cat <- year_cats[y]
              taxonkey <- taxonkeys[t]

              data_thin <- data_thin %>%
                dplyr::filter(group != g) %>%
                dplyr::add_row(year_cat = year_cat,
                               acceptedTaxonKey = taxonkey,
                               n_obs = sum_n_obs,
                               coordinateUncertaintyInMeters = cutoff,
                               decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                               decimalLongitude = data_thin_group$mean_decimalLongitude[1],
                               group = g,
                               geometry = data_thin_group$geometry[1],
                               acceptedScientificName = data_thin_group$acceptedScientificName[1],
                               mean_decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                               mean_decimalLongitude = data_thin_group$mean_decimalLongitude[1])

            }
            if(fun == "n"){
              n_n_obs <- nrow(data_thin_group)

              year_cat <- year_cats[y]
              taxonkey <- taxonkeys[t]

              data_thin <- data_thin %>%
                dplyr::filter(group != g) %>%
                dplyr::add_row(year_cat = year_cat,
                               acceptedTaxonKey = taxonkey,
                               n_obs = n_n_obs,
                               coordinateUncertaintyInMeters = cutoff,
                               decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                               decimalLongitude = data_thin_group$mean_decimalLongitude[1],
                               group = g,
                               geometry = data_thin_group$geometry[1],
                               acceptedScientificName = data_thin_group$acceptedScientificName[1],
                               mean_decimalLatitude = data_thin_group$mean_decimalLatitude[1],
                               mean_decimalLongitude = data_thin_group$mean_decimalLongitude[1])
            }
          }
        }
      }




      # rbind data for each year_cat
      if(y == 1){
        data_thin_all_year <- data_thin
      } else {
        data_thin_all_year <- rbind(data_thin_all_year, data_thin)
      }
    }
    # add back acceptedScientificName
    # get dataframe with acceptedTaxonkey & acceptedScientificNames
    taxa_names <- taxa_data %>%
      sf::st_drop_geometry() %>%
      dplyr::select(acceptedTaxonKey, acceptedScientificName) %>%
      dplyr::distinct() %>%
      dplyr::group_by(acceptedTaxonKey) %>%
      dplyr::add_tally()

    if(nrow(taxa_names) > 1){
      warning(paste0("The data for ", taxonkeys[t], " has more than one acceptedScientificName. >> ", taxa_names$acceptedScientificName[1], " << was used."))
    }

    taxa_names <- taxa_names[1,] %>%
      dplyr::select(-n)

    data_thin_all_year <- data_thin_all_year %>%
      dplyr::select(-acceptedScientificName) %>%
      dplyr::left_join(taxa_names, by = "acceptedTaxonKey")

    # remove group column and mean_decimalLatitude, mean_decimalLongitude
    data_thin_all_year <- data_thin_all_year %>%
      dplyr::select(-group, -mean_decimalLatitude, -mean_decimalLongitude)

    # rbind data for each taxonkey
    if(t == 1){
      data_thin_all <- data_thin_all_year
    } else {
      data_thin_all <- rbind(data_thin_all, data_thin_all_year)
    }
  }
  # return data
  # Convert the data frame to an sf object
  data_thin_all <- sf::st_as_sf(data_thin_all,
                                coords = c("decimalLongitude", "decimalLatitude"),
                                crs = 4326,
                                remove = FALSE)

  return(data_thin_all)
}



