
#' Title
#'
#' @param taxon_key
#' @param basis_of_record
#' @param coord_unc
#' @param identification_verification_status
#'
#' @return
#' @export
#'
#' @examples
#' taxon_key <- c(2865504, 5274858, 527485800)
#' taxon_key <- c(2865504)
#' taxon_key<-4980984 #no occurrences in gbif
#' taxon_key<-4377179 #genus
#' zip_file <- "./<path to zip_file>/0001221-210914110416597.zip"

get_gbif_data <- function(taxon_key,
                          basis_of_record,
                          coord_unc,
                          identification_verification_status
) {

  #-----------------------------------------
  # 1. Test function arguments
  #-----------------------------------------

  # Test that taxon_key is provided
  assertthat::assert_that(!is.null(taxon_key),
                          msg = paste(
                            "taxon_key is missing."
                          )
  )

  # Test that taxon_key is of class numeric
  if (!is.null(taxon_key)) {
    assertthat::assert_that(is.numeric(taxon_key),
                            msg = "taxon_key should be of class numeric."
    )
  }

  # Test that taxon_key (when provided) matches GBIF backbone and that it is an accepted taxon key
  if (!is.null(taxon_key)) {
    taxon_key <- as.integer(taxon_key)
    taxon_key_df <- as.data.frame(taxon_key)

    mapped_taxa <- purrr::map_dfr(
      taxon_key_df$taxon_key,
      ~ {
        tryCatch(
          {
            data <- rgbif::name_usage(key = .x)$data
            if (length(data) == 0) {
              stop("No match with the GBIF backbone found")
            }
            data
          },
          error = function(e) {
            NULL
          }
        )
      }
    )

    assertthat::assert_that(nrow(mapped_taxa)==length(taxon_key),
                            msg=paste0("The following taxon key(s) could not be matched with the GBIF backbone: ",taxon_key[!taxon_key %in% mapped_taxa$key])
    )

    not_accepted <- mapped_taxa %>%
      dplyr::filter(taxonomicStatus !="ACCEPTED")

    if (nrow(not_accepted)!=0) {
      warning(paste0("The following taxon key(s) are not considered accepted taxa in the GBIF backbone: "
                     ,taxon_key[!taxon_key %in% not_accepted$key], ".")
              )
    }
    remove(mapped_taxa)
    remove(not_accepted)
    remove(taxon_key_df)

  }




  #-----------------------------------------
  #           2. Download data
  #-----------------------------------------

  gbif_download <- rgbif::occ_download(
    rgbif::pred_in("taxonKey", taxon_key),
    rgbif::pred("hasCoordinate", TRUE),
    rgbif::pred_gt("year", 1900),
    user = rstudioapi::askForPassword("GBIF username"),
    pwd = rstudioapi::askForPassword("GBIF password"),
    email = rstudioapi::askForPassword("Email address for notification"))

  #Follow the status of the download
  rgbif::occ_download_wait(gbif_download)

  #Retrieve downloaded records
  gbif_data <- rgbif::occ_download_get(gbif_download,overwrite = TRUE) %>%
    rgbif::occ_download_import()

  #Retrieve citation of downloaded dataset
  print(rgbif::gbif_citation(rgbif::occ_download_meta(gbif_download))$download)

  #Check that download provided data
  assertthat::assert_that(nrow(gbif_data)!=0,
                          msg = paste("No occurrence records for taxonkey(s) ",
                            paste(taxon_key, collapse = ", "),
                            " were found on GBIF.")
  )

  #-----------------------------------------
  #         3  . clean data
  #-----------------------------------------
  gbif_data <- gbif_data %>%
    dplyr::select("acceptedTaxonKey",
                  "species",
                  "decimalLatitude",
                  "decimalLongitude",
                  "establishmentMeans",
                  "coordinateUncertaintyInMeters",
                  "basisOfRecord",
                  "taxonRank",
                  "taxonomicStatus",
                  "genus",
                  "specificEpithet",
                  "eventDate",
                  "occurrenceStatus",
                  "gbifID",
                  "year",
                  "countryCode",
                  "identificationVerificationStatus")

  #If the coord_unc argument was (not) provided:
  if(is.null(coord_unc)){
    coord_unc <- max(gbif_data$coordinateUncertaintyInMeters, na.rm = TRUE)
  }

  #If the coord_unc argument was (not) provided:
  if(is.null(identification_verification_status)){
    identificationVerificationStatus_to_discard <- c(
      "unverified",
      "unvalidated",
      "not validated",
      "under validation",
      "not able to validate",
      "control could not be conclusive due to insufficient knowledge",
      "uncertain",
      "unconfirmed",
      "unconfirmed - not reviewed",
      "validation requested"
    )
  }

  # If the basis_of_record argument was (not) provided:
  if(is.null(basis_of_record)){
    basis_of_record <- c(
      "OBSERVATION",
      "HUMAN_OBSERVATION",
      "MATERIAL_SAMPLE",
      "LITERATURE",
      "PRESERVED_SPECIMEN",
      "UNKNOWN",
      "MACHINE_OBSERVATION"
    )

    warning(paste0(
      "The basis_of_record parameter was not provided, the following",
      " record types will be included: ", paste(basis_of_record, collapse = ", ")
    ))
  }

# Create a dataframe with the accepted species name and corresponding acceptedTaxonkey
# This part is only useful when the species column can be different from the name of the accepted species
  species_info <- gbif_data %>%
    dplyr::filter(.data$taxonRank %in% c("SPECIES", "VARIETY"),
                  .data$taxonomicStatus == "ACCEPTED") %>%
    dplyr::distinct(.data$acceptedTaxonKey,
                    .data$genus,
                    .data$specificEpithet) %>%
    dplyr::mutate(acceptedScientificName_2 = paste(.data$genus, .data$specificEpithet)) %>%
    dplyr::rename(acceptedTaxonKey_2 = "acceptedTaxonKey") %>%
    dplyr::distinct(.data$acceptedTaxonKey_2, .data$acceptedScientificName_2)


  data_cleaned <- gbif_data %>%
    dplyr::left_join(species_info, by = c("species" = "acceptedScientificName_2")) %>%
    dplyr::mutate(acceptedTaxonKey = .data$acceptedTaxonKey_2) %>%
    dplyr::filter(
      !is.na(.data$acceptedTaxonKey), #Remove occurrences of species for which the species column does not correspond accepted species (ASN_2)
      !is.na(.data$eventDate),
      !is.na(.data$decimalLatitude),
      .data$eventDate >= "1901-01-01",
      .data$basisOfRecord %in% basis_of_record,
      .data$coordinateUncertaintyInMeters <= coord_unc |
        is.na(.data$coordinateUncertaintyInMeters),
      .data$occurrenceStatus == "PRESENT") %>%
    dplyr::select("gbifID",
                  "year",
                  "acceptedTaxonKey",
                  "decimalLatitude",
                  "decimalLongitude",
                  "coordinateUncertaintyInMeters",
                  "countryCode") %>%
    dplyr::mutate(year_cat = dplyr::case_when(year <= 1925 ~ "1901-1925",
                                              year <= 1950 ~ "1926-1950",
                                              year <= 1975 ~ "1951-1975",
                                              year <= 2000 ~ "1976-2000",
                                              year > 2000 ~ "2001-2025",
                                              TRUE ~ NA_character_)) %>%
    dplyr::group_by(.data$year_cat,
                    .data$acceptedTaxonKey,
                    .data$decimalLatitude,
                    .data$decimalLongitude) %>%
    dplyr::summarize(n_obs = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(species_info, by = c("acceptedTaxonKey" = "acceptedTaxonKey_2")) %>%
    dplyr::rename(acceptedScientificName = "acceptedScientificName_2")

  assertthat::assert_that(nrow(data_cleaned)!=0,
      msg=paste0("No useful data for ",
             paste(taxon_key, collapse = ","),
             paste(" left after filtering.",
                   "Try omiting or changing the filter settings.")
      )
    )


  remove(gbif_data)

  data_sf <- data_cleaned %>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)
}