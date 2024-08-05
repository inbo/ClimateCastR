data_stitching <- function(x,
                           y,
                           scientific_name
) {

  #-----------------------------------------
  # 1. Test function arguments
  #-----------------------------------------

  # Test that x is of class data.frame
  assertthat::assert_that(inherits(x, "data.frame"),
                          msg = paste( "x is not of class data.frame."
                          )
  )

  #Test that y is of class data.frame
  assertthat::assert_that(inherits(y, "data.frame"),
                          msg = paste( "y is not of class data.frame."
                          )
  )

  #Test that scientific name is of class character
  assertthat::assert_that(is.character(scientific_name),
                          msg = paste( "scientific_name is not of class character."
                          )
  )


