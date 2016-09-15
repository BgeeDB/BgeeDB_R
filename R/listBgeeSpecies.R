#' @title List species in the Bgee database and the available data types for each of them
#'
#' @param release A character specifying a targeted release number. In the form "Release.subrelease" or "Release_subrelease", e.g., "13.2" or 13_2". If not specified, the latest release is used.
#'
#' @return A data frame with species Id, genus name, species name, common name and data type availability for targeted Bgee release
#'
#' @examples{
#'  listBgeeSpecies()
#' }
#'
#' @author Julien Roux
#' @export

listBgeeSpecies <- function(release=NULL){
  allReleases <- .getRelease()
  if (length(release)==0) {
    release <- gsub("\\.", "_", allReleases$release[1])
  } else if (length(release)==1){
    # In case the release number is written with a dot
    release <- gsub("\\.", "_", release)
    # test if required release exists
    if (sum(allReleases$release == gsub("_", ".", release))!=1){
      stop("ERROR: The specified release number is invalid")
    }
  } else {
    stop("ERROR: The specified release number is invalid.")
  }

  # Creating the webservice URL to get all species information
  host <- allReleases$TopAnat.URL[allReleases$release == gsub("_", ".", release)]
  myurl <- paste0(host, "?page=dao&action=org.bgee.model.dao.api.species.SpeciesDAO.getAllSpecies&display_type=tsv&attr_list=id&attr_list=genus&attr_list=species_name&attr_list=common_name")
  # TO DO: probably some parameters to add to this URL to show data types

  ## Set the internet.info to 2 to have less verbose output (only reports critical warnings)
  options(internet.info=2)
  ## Set the timeout option to 600 seconds to let some time to the server to send data (default is 60s)
  options(timeout = 600)

  ## Query webservice
  cat(paste0("Query URL successfully built (", myurl,")\nSubmitting URL to Bgee webservice\n  "))
  download.file(myurl, destfile = file.path(getwd(), "allSpecies.tsv"))

  ## Read 5 last lines of file: should be empty indicating success of data transmission
  ## We cannot use a system call to UNIX command since some user might be on Windows
  tmp <- tail(read.table(file.path(getwd(), "allSpecies.tsv"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE), n=5)
  if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
    ## The file transfer was successful!
    cat(paste0("Query to Bgee webservice successful!\n"))
    allSpecies <- read.table(file.path(getwd(), "allSpecies.tsv"), header=TRUE, sep="\t", blank.lines.skip=TRUE, as.is=TRUE)
    ## Remove temporary file
    file.remove(file.path(getwd(), "allSpecies.tsv"))
  } else {
    ## delete the temporary file
    file.remove(file.path(getwd(), "allSpecies.tsv"))
    stop(paste0("ERROR: The queried file is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
  }
  return(allSpecies)
}
## TO DO: sort species based on speciesId? Genus/species name?
