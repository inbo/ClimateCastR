#' Prepare species occurrence data from GBIF for climate casting
#'
#' @param gbif_data A data frame holding occurrence records downloaded from GBIF.
#' @param basis_of_record Optional character indicating the basisOfRecord types to be included in the data.
#' If NULL, the default, occurrences with the following basisOfRecord will be kept:  "OBSERVATION", "HUMAN_OBSERVATION",
#' "MATERIAL_SAMPLE", "LITERATURE", "PRESERVED_SPECIMEN", "UNKNOWN", and "MACHINE_OBSERVATION".
#' @param coord_unc Optional numeric indicating the maximal coordinate uncertainty (m) an
#' observation can have to be included in the data.
#' If NULL, the default, all occurrences will be kept regardless of their coordinate uncertainty.
#' @param identification_verification_status Optional character or a character vector indicating the identificationVerificationStatus of occurrence records that will be kept.
#' If NULL, the default, all occurrences will be kept except those with the following identificationVerificationStatus: "unverified", "unvalidated", "not validated", "under validation", "not able to validate", "control could not be conclusive due to insufficient knowledge",  "Control could not be conclusive due to insufficient knowledge", "1","uncertain", "unconfirmed", "Douteux", "Invalide", "Non réalisable", "verification needed" , "Probable", "unconfirmed - not reviewed", "validation requested".
#'
#' @return  An sf data frame holding occurrence records that are ready to be used for climate casting.
#' @export
#'
#' @examples
data_prep<-function(gbif_data,
                    basis_of_record=NULL,
                    coord_unc=NULL,
                    identification_verification_status=NULL){


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
                  "identificationVerificationStatus",
                  "issue")

  #If the coord_unc argument was not provided:
  if(is.null(coord_unc)){
    coord_unc <- max(gbif_data$coordinateUncertaintyInMeters, na.rm = TRUE)
  }

  #If the identification_verification_status was not provided:
  if(is.null(identification_verification_status)){

    identificationverificationstatus_to_discard<- c(
      "unverified",
      "unvalidated",
      "not validated",
      "under validation",
      "not able to validate",
      "control could not be conclusive due to insufficient knowledge",
      "Control could not be conclusive due to insufficient knowledge",
      "1",
      "uncertain",
      "unconfirmed",
      "Douteux",
      "Invalide",
      "Non r\u00E9alisable",
      "verification needed" ,
      "Probable",
      "unconfirmed - not reviewed",
      "validation requested"
    )

    identification_verification_status<- gbif_data%>%
      dplyr::filter(!.data$identificationVerificationStatus %in% identificationverificationstatus_to_discard)%>%
      dplyr::distinct(.data$identificationVerificationStatus)%>%
      dplyr::pull()

  }


  # If the basis_of_record argument was not provided:
  # We don't take into account `FOSSIL SPECIMEN` and `LIVING SPECIMEN`, which can have misleading location information (see be alientaxa cube)
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

  #First filtering step remove data with default geospatial issues, that are not PRESENT, and that don't have coordinates
  #This was done during the download of the get_data_gbif function but has to be specified separately here
  issues_to_discard <- c(
    "ZERO_COORDINATE",
    "COORDINATE_OUT_OF_RANGE",
    "COORDINATE_INVALID",
    "COUNTRY_COORDINATE_MISMATCH"
  )

  gbif_data<-gbif_data %>%
    dplyr::filter(
      !stringr::str_detect(.data$issue,  paste(issues_to_discard, collapse = "|")),
      .data$occurrenceStatus == "PRESENT",
      !is.na(.data$decimalLongitude) & !is.na(.data$decimalLatitude)
    )



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

  #Do filtering
  suppressMessages(
    data_cleaned <- gbif_data %>%
      dplyr::left_join(species_info, by = c("species" = "acceptedScientificName_2")) %>%
      dplyr::mutate(acceptedTaxonKey = .data$acceptedTaxonKey_2) %>%
      CoordinateCleaner::cc_cen(buffer = 2000) %>% # remove country centroids within 2k
      CoordinateCleaner::cc_cap(buffer = 2000) %>% # remove capitals centroids within 2k
      CoordinateCleaner::cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2k
      dplyr::filter(
        !is.na(.data$acceptedTaxonKey), #Remove occurrences of species for which the species column does not correspond accepted species (ASN_2)
        !is.na(.data$eventDate),
        !is.na(.data$decimalLatitude),
        .data$eventDate >= "1901-01-01",
        .data$basisOfRecord %in% basis_of_record,
        .data$identificationVerificationStatus %in% identification_verification_status,
        .data$coordinateUncertaintyInMeters <= coord_unc |
          is.na(.data$coordinateUncertaintyInMeters)) %>%
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
  )

  assertthat::assert_that(nrow(data_cleaned)!=0,
                          msg=paste0("No useful data left after filtering.Try omiting or changing the filter settings.")
  )


  #Remove occurrences in the ocean, and check that these are not the majority of records
  #NOTE, this still leaves quite some occurrence records very close to the coast!
  suppressMessages(suppressWarnings(
    data_cleaned_no_sea<-data_cleaned%>%
      CoordinateCleaner::cc_sea(verbose=FALSE)
  )
  )

  assertthat::assert_that(nrow(data_cleaned_no_sea) > (nrow(data_cleaned)/2),
                          msg="The majority of occurrence records fall within the ocean. Climate matching is not possible for marine species."
  )

  #-------------------------------------------------
  #      4. Save as an sf dataframe and return
  #-------------------------------------------------
  data_sf <- data_cleaned_no_sea%>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)

  return(data_sf)

}


