#' @title List Bgee releases available to use with BgeeDB package
#'
#' @description Returns information on available Bgee releases, the access URL for FTP and webservice, and the date of release
#'
#' @param release A character specifying a targeted release number. In the form "Release.subrelease" or "Release_subrelease", e.g., "13.2" or 13_2". If not specified, all available releases are shown.
#'
#' @return A data frame with information on Bgee releases.
#'
#' @examples{
#'  listRelease()
#' }
#'
#' @author Julien Roux
#' @export

# Function displaying the user a data frame describing all releases available for Bgee
listRelease <- function(release=NULL){
  allReleases <- .getRelease()
  if (length(release)==1){
    if (sum(allReleases$release == gsub("_", ".", release))==1){
      allReleases <- allReleases[allReleases$release == gsub("_", ".", release), ]
    } else {
      stop("ERROR: The specified release number is invalid.")
    }
  }
  return(allReleases)
}
## TO DO: add some messages to the user

# Function returning a data frame describing all releases available for Bgee
.getRelease <- function(){
  ## query FTP to get file describing all releases
  releaseurl <- 'ftp://ftp.bgee.org/release.tsv'
  try(download.file(releaseurl,
                    destfile=file.path(getwd(), 'release.tsv')
  ),
  silent=FALSE)
  if (!file.exists("release.tsv")){
    stop("ERROR: File describing releases could not be downloaded from FTP.")
  }
  allReleases <- read.table("release.tsv", h=T, sep="\t")
  file.remove(file.path(getwd(), 'release.tsv'))
  return(allReleases)
}
