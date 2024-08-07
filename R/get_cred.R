#' Interactively get & store credentials in the system environment
#'
#' @param x a character with the name of the system environment variable to get
#'   and store if missing
#' @param force a boolean indicating whether to force the user to input the
#' credential
#'
#' @return a character vector containing the value of the system variable
#' @export
#' @importFrom svDialogs dlgInput
#'
#' @examples
#' \dontrun{
#' get_cred("gbif_user")
#' }
get_cred <- function(x, force = FALSE){

  cred <- Sys.getenv(x)

  if(cred == "" | force == TRUE){
    input <- dlgInput(paste0("What is your ", x, "?"))
    cred <- input$res
    do.call(Sys.setenv, as.list(purrr::set_names(cred, x)))
  }
  return(cred)
}
