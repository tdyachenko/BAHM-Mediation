#' Download an instance of bahm to your machine
#'
#' @param path The directory to download bahm to
#' @export
#' @importFrom git2r clone
#' @examples
#' \dontrun{
#'     download_bahm()
#' }
download_bahm <- function(path = getwd()) {
  message("Downloading bahm to ", path, "/bahm")
  path <- file.path(path, "bahm")
  clone("https://github.com/tdyachenko/BAHM-Mediation", local_path = path, branch = "main")
}

#' Run an instance of bahm on your machine
#'
#' @param path The directory to run bahm from
#' @export
#' @importFrom shiny runApp
#' @examples
#' \dontrun{
#'     run_bahm()
#' }
run_bahm <- function(path = getwd()) {
  path <- file.path(path, "bahm", "BAHM")

  if (!file.exists(path)) {
    download_bahm()
  }

  runApp(path)
}
