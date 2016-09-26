#' @title List species in the Bgee database and the available data types for each of them
#'
#' @param release A character specifying a targeted release number. In the form "Release.subrelease" or "Release_subrelease", e.g., "13.2" or 13_2". If not specified, the latest release is used.
#'
#' @param allReleases A data frame with information on all releases. Avoid redownloading this information if .getRelease() already called.
#'
#' @param ordering A numeric indicating which column should be used to sort the data frame. Default NULL, returning unsorted data frame.
#'
#' @return A data frame with species Id, genus name, species name, common name and data type availability for targeted Bgee release
#'
#' @examples{
#'  listBgeeSpecies()
#'  # species present in a specific Bgee release
#'  listBgeeSpecies(release = "13.2")
#'  # in order to order species according to theit taxonomical IDs
#'  listBgeeSpecies(ordering = 1)
#' }
#'
#' @author Julien Roux
#' @export

listBgeeSpecies <- function(release=NULL, ordering=NULL, allReleases=NULL){
  if (length(allReleases)==0) {
    cat("Querying Bgee to get release information...\n")
    allReleases <- .getRelease()
  }
  if (length(release)==0) {
    ## Take latest release
    release <- gsub("\\.", "_", allReleases$release[1])
  } else if (length(release)==1){
    # In case the release number is written with a dot
    release <- gsub("\\.", "_", release)
    # test if required release exists
    if (sum(allReleases$release == gsub("_", ".", release))!=1){
      stop("ERROR: The specified release number is invalid, or is not available for BgeeDB.")
    }
  } else {
    stop("ERROR: The specified release number is invalid.")
  }

  cat(paste0("Building URL to query species in Bgee release ", release, "...\n"))
  host <- allReleases$TopAnat.URL[allReleases$release == gsub("_", ".", release)]
  myurl <- paste0(host, "?page=species&display_type=tsv")

  ## Set the internet.info to 2 to have less verbose output (only reports critical warnings)
  options(internet.info=2)
  ## Set the timeout option to 600 seconds to let some time to the server to send data (default is 60s)
  options(timeout = 600)

  ## Query webservice
  cat(paste0("Submitting URL to Bgee webservice... (", myurl,")\n"))
  download.file(myurl, destfile = file.path(getwd(), "allSpecies.tsv"))

  ## Read 5 last lines of file: should be empty indicating success of data transmission
  ## We cannot use a system call to UNIX command since some user might be on Windows
  tmp <- tail(read.table(file.path(getwd(), "allSpecies.tsv"),
                         header=TRUE,
                         sep="\t",
                         comment.char="",
                         blank.lines.skip=FALSE,
                         as.is=TRUE),
              n=5)
  if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
    ## The file transfer was successful!
    cat(paste0("Query to Bgee webservice successful!\n"))
    allSpecies <- read.table(file.path(getwd(), "allSpecies.tsv"),
                             header=TRUE,
                             sep="\t",
                             blank.lines.skip=TRUE,
                             as.is=TRUE)
    ## Remove temporary file
    file.remove(file.path(getwd(), "allSpecies.tsv"))
  } else {
    ## delete the temporary file
    file.remove(file.path(getwd(), "allSpecies.tsv"))
    stop(paste0("ERROR: The queried file is truncated, ",
                "there may be a temporary problem with the Bgee webservice."))
  }

  if (length(ordering) == 0){
    return(allSpecies)
  } else if(length(ordering) == 1 & is.numeric(ordering) & ordering <= length(allSpecies[1,])){
    return(allSpecies[order(allSpecies[,ordering]),])
  } else {
    cat(paste0("Invalid ordering parameter, returning the data frame unsorted.\n"))
    return(allSpecies)
  }
}
