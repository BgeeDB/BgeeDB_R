# internal function managing download of files. Different approach is used depending on the OS.
# In Windows wininet is used by default but does not allow download longer than 30 seconds.
# In order to remove this timeout problem the package check if libcurl is available. If available
# curl package is used. Otherwise wininet is used and can potentially create timeout.
# in Mac and Linux libcurl is used.
#' @import curl
#' @noMd
#' @noRd

bgee_download_file <- function(url, destfile, quiet = FALSE, mode = NULL) {
  if(.Platform$OS.type == "windows") {
    if(capabilities("libcurl")) {
      if(is.null(mode)) {
        success <- curl_download(url = url, destfile = destfile, quiet = quiet)
      } else{
        success <- curl_download(url = url, destfile = destfile, quiet = quiet, mode = mode)
      }
      if(gsub("\\", "/", success, fixed = TRUE) == destfile | success == destfile) {
        return(0)
      } else {
        return(success)
      }
    } else {
      if(is.null(mode)) {
        warning("Windows OS without compatibility with libcurl. Could result to problems to download files.")
        warning("Please read the installation section of the vignette to solve any potential error.")
        return(download.file(url = url, destfile = destfile, quiet = quiet))
      } else{
        return(download.file(url = url, destfile = destfile, quiet = quiet, mode = mode))
      }
    }
  } else{
    if(is.null(mode)) {
      return(download.file(url = url, destfile = destfile, quiet = quiet))
    } else{
      return(download.file(url = url, destfile = destfile, quiet = quiet, mode = mode))
    }
  }
}