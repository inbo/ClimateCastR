data_stitching <- function(x,
                           y,
                           scientific_name
) {

  #-----------------------------------------
  # 1. Test function arguments
  #-----------------------------------------

  # Test that x is of class sf
  assertthat::assert_that(inherits(x, "sf"),
                          msg = paste( "x is not an sf dataframe."
                          )
  )

  #Test that y is of class sf
  assertthat::assert_that(inherits(y, "sf"),
                          msg = paste( "y is not an sf dataframe."
                          )
  )

  #Test that scientific name is of class character
  assertthat::assert_that(is.character(scientific_name),
                          msg = paste( "scientific_name is not of class character."
                          )
  )

  #Test that all necessary columns are present
  columnnames<-c("acceptedTaxonKey",
                 "decimalLatitude",
                 "decimalLongitude",
                 "coordinateUncertaintyInMeters",
                 "acceptedScientificName",
                 "year_cat",
                 "n_obs",
                 "geometry")

  missing_in_x<- setdiff(names(x),columnnames)

  missing_in_y<- setdiff(names(y),columnnames)

  assertthat::assert_that(columnnames %in% names(x),
                          msg = paste( "The following columns are missing in", x,":", missing_in_x
                          )
  )

  assertthat::assert_that(columnnames %in% names(y),
                          msg = paste( "The following columns are lacking in", y,":", missing_in_y
                          )
  )


  #-----------------------------------------
  # 2. Data prep
  #-----------------------------------------
  x$acceptedScientificName<-scientific_name
  y$acceptedScientificName<-scientific_name

  joined_data<-rbind(x,y)

  return(joined_data)
}
