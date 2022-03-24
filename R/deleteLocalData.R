#' @title Delete local data for the species of the reference class Bgee object
#'
#' @description This function delete data present in the local database. It can delete data 
#' rather from one datatype or all datatypes of one species.
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species
#' 
#' @param allDataTypes A boolean defining rather data of all datatypes from the selected species
#'  should be deleted. If TRUE, the local sqlite database file is deleted. If FALSE, data
#'  corresponding to the datatype of the `myBgeeObject` object will be deleted from the local
#'  sqlite database. 
#' 
#' @author Julien Wollbrett
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'   data <- getData(bgee, experimentId = "SRP007359")
#'   deleteLocalData(bgee, allDataTypes = TRUE)
#' }
#'
#' @export
#' 
deleteLocalData <- function(myBgeeObject, allDataTypes = FALSE) {
  database_path <- sqlitePath(myBgeeObject = myBgeeObject)
  if(allDataTypes) {
    file.remove(database_path)
    message("Local database deleted successfully for ",
      myBgeeObject$speciesName, " [speciesId:", myBgeeObject$speciesId, "]")
  } else {
    conn <- dbConnect(RSQLite::SQLite(), database_path)
    on.exit(dbDisconnect(conn))
    dbRemoveTable(conn, myBgeeObject$dataType)
    message(myBgeeObject$dataType, " data have been removed from the local sqlite ",
      " database for the species ", myBgeeObject$speciesName, " [speciesId:", 
      myBgeeObject$speciesId, "]")
  }
}

#' @title Delete .rds data of one species coming from an old version of the BgeeDB R package.
#'
#' @description Since Bioconductor 3.11 BgeeDB store data in a local SQLite database allowing
#' to optimize memory usage.
#' This function allows to delete .rds files used to store local data from BgeeDB versions older
#' than Bioconductor 3.11, for the species selected in the reference class Bgee object.
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species
#' 
#' @author Julien Wollbrett
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'   deleteOldData(bgee)
#' }
#'
#' @export
#' 
deleteOldData <- function(myBgeeObject) {
  files <- list.files(path = myBgeeObject$pathToData, pattern = "*.rds", full.names = TRUE)
  file.remove(files)
}