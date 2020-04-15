#' @title Delete local database for the species of the reference class Bgee object
#'
#' @description This function delete all data present in the local database of one species.
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species
#' 
#' @author Julien Wollbrett
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'   deleteLocalData(bgee)
#' }
#'
#' @export
#' 
deleteLocalData <- function(myBgeeObject) {
  database_path <- sqlitePath(myBgeeObject = myBgeeObject)
  file.remove(database_path)
  file.create(database_path)
  message("Local database deleted successfully for ", myBgeeObject$speciesName, " [speciesId:", myBgeeObject$speciesId, "]")
}

#' @title Delete .rds data of one species coming from an old version of the BgeeDB R package.
#'
#' @description Since Bioconductor 3.11 BgeeDB store local data in an SQLite database. This allows to only load wanted data in memory.
#' This function allows to delete .rds files used to store local data in previous versions of the package and not used anymore. This 
#' function delete all .rds files for the species selected in the reference class Bgee object.
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