#' @title List Bgee releases available to use with BgeeDB package
#'
#' @description Returns information on available Bgee releases, the access URL for FTP and webservice, and the date of release
#'
#' @param release A character specifying a targeted release number. In the form "Release.subrelease" or "Release_subrelease", e.g., "13.2" or 13_2". If not specified, all available releases are shown.
#'
#' @return A data frame with information on Bgee releases.
#'
#' @examples{
#'  listBgeeRelease()
#' }
#'
#' @author Julien Roux, Julien Wollbrett
#' @export

# Function displaying the user a data frame describing all releases available for Bgee
listBgeeRelease <- function(release=NULL){
  cat("Downloading release information from Bgee...\n")
  allReleases <- .getBgeeRelease()
  if (length(release)==1){
    if (sum(allReleases$release == gsub("_", ".", release))==1){
      cat(paste0("Only displaying information from targeted release ",
                 gsub("_", ".", release), "\n"))
      allReleases <- allReleases[allReleases$release == gsub("_", ".", release), ]
    } else {
      stop("ERROR: The specified release number is invalid or is not available for BgeeDB.")
    }
  }
  ## Only return the columns of interest to the user
  return(allReleases[, c("release","release.date", "FTP.URL", "TopAnat.URL")])
}

# Function returning a data frame describing all releases available for Bgee
.getBgeeRelease <- function(removeFile=TRUE){
  ## query FTP to get file describing all releases
  #TODO allow download of release from the release_v2 file stored on the ftp
  releaseUrl <- 'ftp://ftp.bgee.org/release_v2.tsv'
  success <- try(download.file(releaseUrl, quiet = TRUE,
                           destfile=file.path(getwd(), 'release.tsv.tmp')))
  if (success == 0 & file.exists(file.path(getwd(), 'release.tsv.tmp'))){
    file.rename(from=file.path(getwd(), 'release.tsv.tmp'),
                to=file.path(getwd(), 'release.tsv'))
    allReleases <- read.table("release.tsv", header=TRUE, sep="\t")
    if (removeFile == TRUE){
      file.remove(file.path(getwd(), 'release.tsv'))
    }
  } else {
    file.remove(file.path(getwd(), 'release.tsv.tmp'))
    stop("ERROR: File describing releases could not be downloaded from FTP.")
  }
  #allReleases <- read.table("release.tsv", header=TRUE, sep="\t")
  sapply(as.character(allReleases$minimumVersionBgeeDB),compareVersion,as.character(packageVersion("BgeeDB")))
  allAvailableReleases <- allReleases[  sapply(as.character(allReleases$minimumVersionBgeeDB),compareVersion,as.character(packageVersion("BgeeDB"))) >= 0,]
  return(allAvailableReleases)
}
