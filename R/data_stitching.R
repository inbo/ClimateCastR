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


