#' @title Retrieve Bgee processed expression values at cell level for single cell data in Bgee.
#'
#' @description This function loads the processed gene count matrix of single cell experiments in Bgee.
#' These data are available for all droplet based single-cell RNA-seq and full length single-cell RNA-seq.
#' These files are stored in the H5AD format.
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species 
#' and data type.
#'
#' @param experimentId Filter allowing to specify one ArrayExpress or GEO accession, e.g., 
#' GSE43721. Can not be null
#' 
#' @param package Name of the R package used to load the H5AD cell data. It can either be \"anndata\"
#' or \"zellkonverter\". If zellkonverter is used then its default reader (R) is used.
#' 
#' @author Julien Wollbrett
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Gallus_gallus", dataType = "sc_droplet_based")
#'   cellProcessedData <- getCellProcessedData(bgee, experimentId = "ERP132579")
#'   
#' }
#'
#' @import zellkonverter anndata HDF5Array
#' @export
#' 
getCellProcessedData <- function(myBgeeObject, experimentId, package = "zellkonverter") {
  if (length(experimentId) != 1) {
    stop("One experimentId has to be provided.")
  }
  if(length(package) != 1 | ! package %in% c("anndata", "zellkonverter")) {
    stop("the R package used to load H5AD cell data should either be anndata or zellkonverter.")
  }
  if (length(myBgeeObject$dataType) != 1 | ! myBgeeObject$dataType %in% c("sc_droplet_based", "sc_full_length")) {
    stop("Cell raw data are only available for single cell data. Please choose \"sc_droplet_based\" or \"sc_full_length\"",
    " as a datatpe of the Bgee object to be able to download such data")
  }
  experimentAnnotation <- getAnnotation(myBgeeObject)$experiment.annotation
  if (!experimentId %in% experimentAnnotation$Experiment.ID) {
    stop("the experiment ID:", experimentId, " provided does not exist in Bgee for the species",
      myBgeeObject$species, ". Please first look at the available annotation using the ",
      "`getAnnotation(myBgeeObject)` function.")
  }
  destFile <- file.path(myBgeeObject$pathToData, paste0(experimentId, "_", myBgeeObject$dataType, ".h5ad"))
  downloadURL <- gsub("EXPIDPATTERN", experimentId, myBgeeObject$experimentH5adUrl)

  if (!file.exists(destFile)) {
    bgee_download_file(url = downloadURL, destfile = destFile, mode = 'wb')
  }
  error_message <- NULL
  experimentH5ad <- NULL
  if (package == "zellkonverter") {
    experimentH5ad <- readH5AD(file = destFile, reader = "R")
  } else if (package == "anndata") {
    tryCatch( { experimentH5ad <- read_h5ad(filename = destFile) },
      error = function(e) {error_message <<- e})
    experimentH5ad <- read_h5ad(filename = destFile)
  }
  if (! is.null(error_message)) {
    stop("an error occured : ", error_message)
  }
  return(experimentH5ad)
}
