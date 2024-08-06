
#Create helper function to clean the data for get_data_gbif, get_data_zip, and get_data_downloadkey
data_clean<-function(gbif_data,
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



  return(data_cleaned_no_sea)

}



#' Download species occurrence data from GBIF and prepare them for climate casting
#'
#' @description
#' This function downloads species occurrence records from GBIF, filters them,
#' and converts them into an sf data frame which can then be used for climate casting.
#'
#' @param taxon_key A numeric value or vector specifying the GBIF taxonKey(s).
#' @param basis_of_record Optional character indicating the basisOfRecord types to be included in the data.
#' If NULL, the default, occurrences with the following basisOfRecord will be kept:  "OBSERVATION", "HUMAN_OBSERVATION",
#' "MATERIAL_SAMPLE", "LITERATURE", "PRESERVED_SPECIMEN", "UNKNOWN", and "MACHINE_OBSERVATION".
#' @param coord_unc Optional numeric indicating the maximal coordinate uncertainty (m) an
#' observation can have to be included in the data.
#' If NULL, the default, all occurrences will be kept regardless of their coordinate uncertainty.
#' @param identification_verification_status Optional character or a character vector indicating the identificationVerificationStatus of occurrence records that will be kept.
#' If NULL, the default, all occurrences will be kept except those with the following identificationVerificationStatus: "unverified", "unvalidated", "not validated", "under validation", "not able to validate", "control could not be conclusive due to insufficient knowledge",  "Control could not be conclusive due to insufficient knowledge", "1","uncertain", "unconfirmed", "Douteux", "Invalide", "Non réalisable", "verification needed" , "Probable", "unconfirmed - not reviewed", "validation requested".
#' @param region_shape Optional character specifying a path to a shape file indicating the region from which occurrence records should be downloaded.
#'
#' @returns An sf data frame holding occurrence records that are ready to be used for climate casting.
#'
#' @export
#'
#' @examples
#'  \dontrun{
#' #provide GBIF taxon_key(s)
#' taxon_key <- c(2865504, 5274858)
#' get_gbif_data(taxon_key)
#' get_gbif_data(taxon_key,
#'               coord_unc = 100,
#'               basis_of_record = "HUMAN_OBSERVATION",
#'               identification_verification_status= c("Probable", "valid","confident","")
#'               )
#'
#' #provide a taxon key and a path to a shapefile
#' taxon_key<-2427091
#' url <- "https://zenodo.org/records/3386224"
#' zipfile <- tempfile(fileext = ".zip")
#' utils::download.file(url, zipfile, mode = "wb")
#' utils::unzip(zipfile, exdir = tempdir())
#' shapefile_path<-file.path(tempdir(), "flanders.shp")
#'
#' get_gbif_data(taxon_key,
#'              region_shape = shapefile_path
#'              )
#'}
#'
#' @author Soria Delva, Sander Devisscher

get_gbif_data <- function(taxon_key,
                          basis_of_record = NULL,
                          coord_unc = NULL,
                          identification_verification_status= NULL,
                          region_shape=NULL
) {

  #-----------------------------------------
  # 1. Test function arguments
  #-----------------------------------------

  # Test that taxon_key is provided
  assertthat::assert_that(!is.null(taxon_key),
                          msg = paste( "taxon_key is missing."
                          )
  )

  # Test that taxon_key is of class numeric
  assertthat::assert_that(is.numeric(taxon_key),
                          msg = "taxon_key should be of class numeric."
  )


  # Test that basis_of_record (when provided) is of class character
  if (!is.null(basis_of_record)) {
    assertthat::assert_that(is.character(basis_of_record),
                            msg = "basis_of_record should be of class character."
    )
  }

  # Test that coord_unc (when provided) is of class numeric
  if (!is.null(coord_unc)) {
    assertthat::assert_that(is.numeric (coord_unc),
                            msg = "coord_unc should be of class numeric."
    )
  }

  # Test that identification_verification_status (when provided) is of class character
  if (!is.null(identification_verification_status)) {
    assertthat::assert_that(is.character(identification_verification_status),
                            msg = "identification_verification_status should be of class character."
    )
  }

  # Test that region_shape (when provided) is of class character
  if (!is.null(region_shape)) {
    assertthat::assert_that(is.character(region_shape),
                            msg = "region_shape should be of class character."
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
      dplyr::filter(.data$taxonomicStatus !="ACCEPTED")

    if (nrow(not_accepted)!=0) {
      warning(paste0("The following taxon key(s) are not considered accepted taxa in the GBIF backbone: "
                     ,taxon_key[!taxon_key %in% not_accepted$key], ".")
      )
    }

  }


  #-----------------------------------------
  #           2. Download data
  #-----------------------------------------

  #Set credentials
  gbif_user <- get_cred("gbif_username")
  gbif_pwd <- get_cred("gbif_password")
  gbif_email <- get_cred("gbif_emailaddress")

  if (!is.null(region_shape)) {
    #read in shapefile
    region_sf<- sf::st_read(region_shape, quiet=TRUE)


    #Place it in the right orientation
    region_sf<- region_sf %>%
      sf::st_transform(crs=4326)%>%
      sf::st_simplify(dTolerance=10) %>% #simplify the polygon because there may be too many points for rgbif
      sf::st_geometry() %>%
      sf::st_as_text()%>%
      wk::wkt() %>%
      wk::wk_orient()

    gbif_download <- rgbif::occ_download(
      rgbif::pred_in("taxonKey", taxon_key),
      rgbif::pred("hasCoordinate", TRUE),
      rgbif::pred_gt("year", 1900),
      rgbif::pred("occurrenceStatus","PRESENT"),
      rgbif::pred("hasGeospatialIssue", FALSE), # remove GBIF default geospatial issues
      rgbif::pred_within(region_sf),
      user = gbif_user,
      pwd = gbif_pwd,
      email = gbif_email)
  } else {
    gbif_download <- rgbif::occ_download(
      rgbif::pred_in("taxonKey", taxon_key),
      rgbif::pred("hasCoordinate", TRUE),
      rgbif::pred_gt("year", 1900),
      rgbif::pred("occurrenceStatus","PRESENT"),
      rgbif::pred("hasGeospatialIssue", FALSE), # remove GBIF default geospatial issues
      user = gbif_user,
      pwd = gbif_pwd,
      email = gbif_email)
  }

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
  #         3. Clean data
  #-----------------------------------------
  gbif_data_clean <- gbif_data %>%
    data_clean(basis_of_record, coord_unc, identification_verification_status)

  #-------------------------------------------------
  #      4. Save as an sf dataframe and return
  #-------------------------------------------------
  data_sf <- gbif_data_clean%>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)

  return(data_sf)
}





#' Import species occurrence data from a zip file and prepare them for climate casting
#'
#' @description
#' This function imports species occurrence data from a zip file, filters them,
#' and converts them into an sf data frame which can then be used for climate casting.
#'
#' @param zip_file Character specifying the path to a zip file holding occurrence records from a previous GBIF download.
#' @param basis_of_record Optional character indicating the basisOfRecord types to be included in the data.
#' If NULL, the default, occurrences with the following basisOfRecord will be kept:  "OBSERVATION", "HUMAN_OBSERVATION",
#' "MATERIAL_SAMPLE", "LITERATURE", "PRESERVED_SPECIMEN", "UNKNOWN", and "MACHINE_OBSERVATION".
#' @param coord_unc Optional numeric indicating the maximal coordinate uncertainty (m) an
#' observation can have to be included in the data.
#' If NULL, the default, all occurrences will be kept regardless of their coordinate uncertainty.
#' @param identification_verification_status Optional character or a character vector indicating the identificationVerificationStatus of occurrence records that will be kept.
#' If NULL, the default, all occurrences will be kept except those with the following identificationVerificationStatus: "unverified", "unvalidated", "not validated", "under validation", "not able to validate", "control could not be conclusive due to insufficient knowledge",  "Control could not be conclusive due to insufficient knowledge", "1","uncertain", "unconfirmed", "Douteux", "Invalide", "Non réalisable", "verification needed" , "Probable", "unconfirmed - not reviewed", "validation requested".
#'
#' @returns An sf data frame holding occurrence records that are ready to be used for climate casting.
#' @export
#'
#' @examples
#' #'  \dontrun{
#' #download zip_file from GBIF
#' zipdir<-tempdir()
#' rgbif::occ_download_get("0001221-210914110416597",
#'                       path = zipdir,
#'                       overwrite = TRUE)
#' zip_file <- file.path(zipdir, "0001221-210914110416597.zip")
#' file.exists(zip_file)
#' get_zip_data(zip_file)
#'}
#'
#' @author Soria Delva, Sander Devisscher

get_zip_data <- function(zip_file,
                         basis_of_record = NULL,
                         coord_unc = NULL,
                         identification_verification_status= NULL
) {

  #-----------------------------------------
  # 1. Test function arguments
  #-----------------------------------------

  # Test that taxon_key is provided
  assertthat::assert_that(!is.null(zip_file),
                          msg = paste(
                            "zip_file is missing."
                          )
  )

  # Test that zip_file is of class character
  assertthat::assert_that(is.character(zip_file),
                          msg = "zip_file should be of class character."
  )



  # Test that the path to the zip file is correct, and the file exists
  assertthat::assert_that(file.exists(zip_file),
                          msg = "The zip file cannot be found. Please check the specified path."
  )

  # Test that basis_of_record (when provided) is of class character
  if (!is.null(basis_of_record)) {
    assertthat::assert_that(is.character(basis_of_record),
                            msg = "basis_of_record should be of class character."
    )
  }

  # Test that coord_unc (when provided) is of class numeric
  if (!is.null(coord_unc)) {
    assertthat::assert_that(is.numeric (coord_unc),
                            msg = "coord_unc should be of class numeric."
    )
  }

  # Test that identification_verification_status (when provided) is of class character
  if (!is.null(identification_verification_status)) {
    assertthat::assert_that(is.character(identification_verification_status),
                            msg = "identification_verification_status should be of class character."
    )
  }


  #-----------------------------------------
  #           2. Download data
  #-----------------------------------------
  expected_columns <- c("acceptedTaxonKey",
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


  zip_data <- readr::read_tsv(unz(zip_file, "occurrence.txt"),
                              col_types = c(decimalLatitude = readr::col_number(),
                                            decimalLongitude = readr::col_number(),
                                            establishmentMeans = readr::col_character(),
                                            acceptedTaxonKey = readr::col_integer()),
                              col_select = expected_columns)

  #Check that all the required columns are present
  assertthat::assert_that(all(expected_columns %in% names(zip_data)),
                          msg = paste("The following columns are missing from the occurrence data:",
                                      paste(expected_columns[!expected_columns %in% names(zip_data)], collapse = ", ")))

  #Check that the zip file contains data
  assertthat::assert_that(nrow(zip_data)!=0,
                          msg = paste("No data were found in the specified zip file.")
  )



  #-----------------------------------------
  #         3  . Clean data
  #-----------------------------------------
  #We will remove data with default geospatial issues, that are not PRESENT, and that don't have coordinates
  #This was done during the download of the get_data_gbif function but has to be specified separately here
  issues_to_discard <- c(
    "ZERO_COORDINATE",
    "COORDINATE_OUT_OF_RANGE",
    "COORDINATE_INVALID",
    "COUNTRY_COORDINATE_MISMATCH"
  )

  zip_data<- zip_data %>%
    dplyr::filter(
      !stringr::str_detect(.data$issue,  paste(issues_to_discard, collapse = "|")),
      .data$occurrenceStatus == "PRESENT",
      !is.na(.data$decimalLongitude) & !is.na(.data$decimalLatitude)
    )


  zip_data_clean <- zip_data %>%
    data_clean(basis_of_record, coord_unc, identification_verification_status)

  #-------------------------------------------------
  #      4. Save as an sf dataframe and return
  #-------------------------------------------------
  zip_data_sf <- zip_data_clean%>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)

  return(zip_data_sf)
}




#' Import species occurrence data from a previous GBIF download and prepare them for climate casting
#'
#' This function imports species occurrence data from a previous GBIF download, filters them,
#' and converts them into an sf data frame which can then be used for climate casting.
#'
#' @param downloadkey Character specifying the download key of a previous GBIF download.
#' @param basis_of_record Optional character indicating the basisOfRecord types to be included in the data.
#' If NULL, the default, occurrences with the following basisOfRecord will be kept:  "OBSERVATION", "HUMAN_OBSERVATION",
#' "MATERIAL_SAMPLE", "LITERATURE", "PRESERVED_SPECIMEN", "UNKNOWN", and "MACHINE_OBSERVATION".
#' @param coord_unc Optional numeric indicating the maximal coordinate uncertainty (m) an
#' observation can have to be included in the data.
#' If NULL, the default, all occurrences will be kept regardless of their coordinate uncertainty.
#' @param identification_verification_status Optional character or a character vector indicating the identificationVerificationStatus of occurrence records that will be kept.
#' If NULL, the default, all occurrences will be kept except those with the following identificationVerificationStatus: "unverified", "unvalidated", "not validated", "under validation", "not able to validate", "control could not be conclusive due to insufficient knowledge",  "Control could not be conclusive due to insufficient knowledge", "1","uncertain", "unconfirmed", "Douteux", "Invalide", "Non réalisable", "verification needed" , "Probable", "unconfirmed - not reviewed", "validation requested".
#'
#' @returns An sf data frame holding occurrence records that are ready to be used for climate casting.
#' @export
#'
#' @author Soria Delva, Sander Devisscher
#' @examples
#' #Provide a downloadkey
#' downloadkey<-"0001221-210914110416597"
#' get_downloadkey_data(downloadkey,
#'                      coord_unc= 1000,
#'                      identification_verification_status= c("","valid"))
#'

get_downloadkey_data <- function(downloadkey,
                                 basis_of_record = NULL,
                                 coord_unc = NULL,
                                 identification_verification_status= NULL
) {

  #-----------------------------------------
  # 1. Test function arguments
  #-----------------------------------------

  # Test that downloadkey is provided
  assertthat::assert_that(!is.null(downloadkey),
                          msg = paste(
                            "downloadkey is missing."
                          )
  )

  # Test that downloadkey is of class character
  assertthat::assert_that(is.character(downloadkey),
                          msg = "downloadkey should be of class character."
  )

  # Test that basis_of_record (when provided) is of class character
  if (!is.null(basis_of_record)) {
    assertthat::assert_that(is.character(basis_of_record),
                            msg = "basis_of_record should be of class character."
    )
  }

  # Test that coord_unc (when provided) is of class numeric
  if (!is.null(coord_unc)) {
    assertthat::assert_that(is.numeric (coord_unc),
                            msg = "coord_unc should be of class numeric."
    )
  }

  # Test that identification_verification_status (when provided) is of class character
  if (!is.null(identification_verification_status)) {
    assertthat::assert_that(is.character(identification_verification_status),
                            msg = "identification_verification_status should be of class character."
    )
  }


  #-----------------------------------------
  #           2. Download data
  #-----------------------------------------

  #Retrieve downloaded records
  downloadkey_data <- rgbif::occ_download_get(downloadkey,overwrite = TRUE) %>%
    rgbif::occ_download_import()

  #Retrieve citation of downloaded dataset
  print(rgbif::gbif_citation(rgbif::occ_download_meta(downloadkey))$download)


  #Check that the zip file contains data
  assertthat::assert_that(nrow(downloadkey_data)!=0,
                          msg = paste("No data were found using the specified downloadkey.")
  )



  #-----------------------------------------
  #         3  . Clean data
  #-----------------------------------------

  #We will remove data with default geospatial issues, that are not PRESENT, and that don't have coordinates
  #This was done during the download of the get_data_gbif function but has to be specified separately here
  issues_to_discard <- c(
    "ZERO_COORDINATE",
    "COORDINATE_OUT_OF_RANGE",
    "COORDINATE_INVALID",
    "COUNTRY_COORDINATE_MISMATCH"
  )

  downloadkey_data<- downloadkey_data %>%
    dplyr::filter(
      !stringr::str_detect(.data$issue,  paste(issues_to_discard, collapse = "|")),
      .data$occurrenceStatus == "PRESENT",
      !is.na(.data$decimalLongitude) & !is.na(.data$decimalLatitude)
    )

    downloadkey_data_clean <- downloadkey_data %>%
    data_clean(basis_of_record, coord_unc, identification_verification_status)

  #-------------------------------------------------
  #      4. Save as an sf dataframe and return
  #-------------------------------------------------
  downloadkey_data_sf <- downloadkey_data_clean%>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)

  return(downloadkey_data_sf)
}


