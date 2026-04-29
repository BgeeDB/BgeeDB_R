#' @title List species in the Bgee database and the available data types for each of them
#'
#' @description Returns information on available species in Bgee
#'
#' @param release A character specifying a targeted release number. In the form "Release.subrelease" or "Release_subrelease", e.g., "13.2" or 13_2". If not specified, the latest release is used.
#'
#' @param allReleases A data frame with information on all releases. Avoid redownloading this information if .getBgeeRelease() already called.
#'
#' @param ordering A numeric indicating the number of the column which should be used to sort the data frame. Default NULL, returning unsorted data frame.
#'
#' @param removeFile Boolean indicating whether the downloaded file should be deleted. Default to TRUE.
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

listBgeeSpecies <- function(release=NULL, ordering=NULL, allReleases=NULL, removeFile=TRUE){
  
  OLD_WEBSERVICE_VERSION = '13.2'
  NEW_VERSION_FROM_FTP_FILE = '15.2'

  if (length(allReleases)==0) {
    cat("\nQuerying Bgee to get release information...\n")
    allReleases <- .getBgeeRelease()
  }
  if (length(release)==0) {
    ## Take latest release
    release <- gsub("\\.", "_", allReleases$release[1])
  } else if (length(release)==1){
    # In case the release number is written with a dot
    release <- gsub("\\.", "_", release)
    # test if required release exists
    if (sum(as.numeric(allReleases$release) == as.numeric(gsub("_", ".", release)))!=1){
      stop("ERROR: The specified release number is invalid, or is not available for BgeeDB.")
    }
  } else {
    stop("ERROR: The specified release number is invalid.")
  }
  
  speciesFilePath <- paste0(getwd(), "/species_Bgee_", release, ".tsv")
  cat(paste0("\nBuilding URL to query species in Bgee release ", release, "...\n"))
  
  myUrl <- allReleases$TopAnat.URL[as.numeric(allReleases$release) == as.numeric(gsub("_", ".", release))]
  
  #create the url depending on selected release of Bgee. It uses a webservice for releases before 15.2 and a file starting from 15.2
  if(compareVersion(gsub("_", ".",release), OLD_WEBSERVICE_VERSION) > 0){
    if (compareVersion(gsub("_", ".",release), NEW_VERSION_FROM_FTP_FILE) >= 0){
      myUrl <- file.path(paste0("https://www.bgee.org/ftp/bgee_v", release), "rPackageSpeciesInfo.tsv")
    } else {
      myUrl <- paste0(myUrl, "?page=r_package&action=get_all_species&display_type=tsv&source=BgeeDB_R_package&source_version=", packageVersion("BgeeDB"))
    }
  } else {
    myUrl <- paste0(myUrl, "?page=species&display_type=tsv&source=BgeeDB_R_package&source_version=", packageVersion("BgeeDB"))
  }

  ## Set the internet.info to 2 to have less verbose output (only reports critical warnings)
  options(internet.info=2)
  ## Set the timeout option to 600 seconds to let some time to the server to send data (default is 60s)
  options(timeout = 600)

  ## Query webservice
  cat(paste0("\ndownloading Bgee species info... (", myUrl,")\n"))
  bgee_download_file(url = myUrl, destfile = speciesFilePath)
  ## Read 5 last lines of file: should be empty indicating success of data transmission
  ## We cannot use a system call to UNIX command since some user might be on Windows
  tmp <- tail(read.table(speciesFilePath,
                         header=TRUE,
                         sep="\t",
                         comment.char="",
                         blank.lines.skip=FALSE,
                         as.is=TRUE),
              n=5)
  if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
    ## The file transfer was successful!
    cat(paste0("\nDownload of species information successful!\n"))
    allSpecies <- read.table(speciesFilePath,
                             header=TRUE,
                             sep="\t",
                             blank.lines.skip=TRUE,
                             as.is=TRUE)
    if (removeFile == TRUE){
      ## Remove temporary file
      file.remove(speciesFilePath)
    }
  } else {
    ## delete the temporary file
    file.remove(speciesFilePath)
    stop(paste0("ERROR: The queried file is truncated, ",
                "there may be a temporary problem with the Bgee webservice."))
  }

  if (length(ordering) == 0){
    return(allSpecies)
  } else if(length(ordering) == 1 & is.numeric(ordering) & ordering <= length(allSpecies[1,])){
    return(allSpecies[order(allSpecies[,ordering]),])
  } else {
    cat(paste0("\nWARNING: Invalid ordering parameter, returning the data frame unsorted.\n"))
    return(allSpecies)
  }
}
