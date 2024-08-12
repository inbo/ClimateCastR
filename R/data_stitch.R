
#' Combine occurrence data from different species into a single species complex and prepare this dataset for climate casting.
#'
#' @param x An sf dataframe, resulting from either get_gbif_data(), get_zip_data(), or get_downloadkey_data().
#' @param to_stitch A numeric vector of two or more taxon keys, present in x under "acceptedTaxonKey", that should be combined into a species complex.
#' @param taxa_complex_name A character value, specifying the name that will be assigned to the species complex.
#' @param taxa_complex_key An optional numeric value, specifying the acceptedTaxonKey that will be assigned to the species complex. If not provided,the taxon key of one of the species in the complex is used.
#'
#' @return A sf data frame, holding occurrence records of a species complex, ready to be used for climate casting.
#' @export
#'
#' @examples
#'  \dontrun{
#' taxon_key <- c(2865504,2891932,2891930)
#' df<-get_gbif_data(taxon_key,
#'                   path = tempdir()
#'                   )
#' df_stitched<- data_stitch(df,
#'                           to_stitch = c(2891932,2891930),
#'                           taxa_complex_name= "Casuarina cunninghamiana/equisetifolia"
#'                           )
#'           }
data_stitch <- function(x,
                        to_stitch,
                        taxa_complex_name,
                        taxa_complex_key = NULL
) {

  #-----------------------------------------
  # 1. Test function arguments
  #-----------------------------------------

  # Test that x is of class sf
  assertthat::assert_that(inherits(x, "sf"),
                          msg = paste( "x is not an sf dataframe."
                          )
  )


  #Test that contained_taxa is of class numeric
  assertthat::assert_that(is.numeric(to_stitch),
                          msg = paste( "to_stitch is not of class numeric."
                          )
  )

  #Test that taxa_complex_name is of class character
  assertthat::assert_that(is.character(taxa_complex_name),
                          msg = paste( "taxa_complex_name is not of class character."
                          )
  )

  #Test that, when provided, taxa_complex_key is of class integer
  if (!is.null(taxa_complex_key)) {
    assertthat::assert_that(is.numeric(taxa_complex_key),
                            msg = "taxa_complex_key should be of class numeric."
    )
  }

  #Test that all necessary columns are present
  columnnames<-c("acceptedTaxonKey",
                 "decimalLatitude",
                 "decimalLongitude",
                 "coordinateUncertaintyInMeters",
                 "acceptedScientificName",
                 "year_cat",
                 "n_obs",
                 "geometry")
  colnames_x <-names(x)

  missing_in_x<- setdiff(columnnames, colnames_x)


  assertthat::assert_that(all(columnnames %in% colnames_x),
                          msg = paste( "The following columns are missing in", deparse(substitute(x)),":", paste(missing_in_x, collapse = ", "))
  )

  #-----------------------------------------
  # 2. Data prep
  #-----------------------------------------

  #Create a dataset with data that should be stitched and one with data that should not be touched
  data_to_keep<- x %>%
    dplyr::filter(!acceptedTaxonKey %in% to_stitch)

  data_to_stitch<- x %>%
    dplyr::filter(acceptedTaxonKey %in% to_stitch)

  # Add a common taxonkey and scientific name to the data_to_stich dataset
  data_to_stitch$acceptedScientificName<-taxa_complex_name

  if(!is.null (taxa_complex_key)) {
  data_to_stitch$acceptedTaxonKey<- taxa_complex_key
  }else{
  data_to_stitch$acceptedTaxonKey<- dplyr::first(unique(data_to_stitch$acceptedTaxonKey))
  }

  #Recalculate the n_obs per year_category for this species complex
  data_to_stitch<-data_to_stitch %>%
    dplyr::group_by(.data$year_cat,
                    .data$acceptedScientificName,
                    .data$acceptedTaxonKey,
                    .data$decimalLatitude,
                    .data$decimalLongitude) %>%
    dplyr::mutate(coordinateUncertaintyInMeters = ifelse(all(is.na(.data$coordinateUncertaintyInMeters)),
                                                         NA_real_,
                                                         min(.data$coordinateUncertaintyInMeters, na.rm = TRUE)))%>%
    dplyr::summarize(n_obs = dplyr::n(),
                     coordinateUncertaintyInMeters = dplyr::first(.data$coordinateUncertaintyInMeters)) %>%
    dplyr::ungroup() %>%
    dplyr::relocate(acceptedScientificName, .after=coordinateUncertaintyInMeters) #reorder columns

  #Combine the two dataframes and return
  stitched_data<-rbind(data_to_keep, data_to_stitch)

  return(stitched_data)
}
