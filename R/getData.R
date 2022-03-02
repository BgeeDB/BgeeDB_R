#' @title Retrieve Bgee RNA-seq or Affymetrix data.
#'
#' @description This function loads the quantitative expression data and presence calls 
#' for samples available from Bgee (rna_seq, affymetrix).
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species 
#' and data type.
#'
#' @param experimentId Filter allowing to specify one or more ArrayExpress or GEO accession, e.g., 
#' GSE43721. Default is NULL: takes all available experiments for targeted species and data type.
#' 
#' @param sampleId Filter allowing to specify one or more sample ID. Depending on the selected 
#' datatype this sample IDs can correspond to Chip IDs (affymetrix) or RNA-Seq library IDs (rna_seq). 
#' Default is NULL: takes all available samples for targeted species and data type.
#'
#' @param anatEntityId Filter allowing to specify one or more anatomical entity IDs from the UBERON 
#' ontology (http://uberon.github.io/). Default is NULL: takes all available anatomical entities for 
#' targeted species and data type.
#'
#' @param stageId Filter allowing to specify one or more developmental stage IDs from Developmental 
#' Stage Ontology (https://github.com/obophenotype/developmental-stage-ontologies). Default is 
#' NULL: takes all available developmental stages for targeted species and data type.
#' 
#' @return Return a dataframe containing all Bgee processed expression data from the selected species 
#' and datatype using specified filters with operator AND.
#'
#' @author Julien Wollbrett, Andrea Komljenovic and Julien Roux.
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'   dataMouse <- getData(bgee)
#'   dataMouseGSE43721 <- getData(bgee, experimentId = "GSE43721")
#'   dataMouseVariousFilters <- getData(bgee, experimentId = c("GSE43721", "GSE36026"), 
#'                              anatEntityId = c("UBERON:0002107", "UBERON:0000956", "UBERON:0002048"))
#' }
#'
#' @import RSQLite
#' @export
#' 
getData <- function(myBgeeObject, experimentId = NULL, sampleId = NULL, anatEntityId = NULL, stageId = NULL) {
  check_object(myBgeeObject)
  import_data(myBgeeObject = myBgeeObject, experimentId = experimentId, sampleId = sampleId, anatEntityId = anatEntityId, stageId = stageId)
  return(query_data(myBgeeObject, experimentId, sampleId, anatEntityId, stageId))
}

check_object = function(myBgeeObject, experimentId = NULL){
  ## check that the Bgee object is valid
  if (length(myBgeeObject$quantitativeData) == 0 ){
    stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object was correctly built.")
  } else {
    ## Check that download of quantitative data is possible for targeted species and data type
    ## Try to output error message that helps a bit the user
    if (myBgeeObject$quantitativeData != TRUE){
      if (length(myBgeeObject$dataType) > 1){
        stop("ERROR: downloading quantitative data is only possible if a single data type (\"rna_seq\" or \"affymetrix\") is specified in the input Bgee class object.")
      } else if (length(myBgeeObject$dataType) == 1 & !(myBgeeObject$dataType %in% c('rna_seq','affymetrix'))){
        stop("ERROR: downloading quantitative data is not possible for the species and data type specified in the input Bgee class object. Maybe the specified data type is not available for the targeted species in Bgee? See listBgeeSpecies() for details on data types availability for each species.")
      } else {
        stop("ERROR: downloading quantitative data is not possible for the species and data type specified in the input Bgee class object.")
      }
    } else if (length(myBgeeObject$experimentUrl) == 0 | length(myBgeeObject$allExperimentsUrl) == 0 | length(myBgeeObject$dataType) == 0 | length(myBgeeObject$pathToData) == 0){
      stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object was correctly built.")
    }
  }
}

## extract the global path to the sqliteDB
sqlitePath <- function(myBgeeObject){
  sqlite_file <- file.path(myBgeeObject$pathToData, paste0(myBgeeObject$speciesName, myBgeeObject$sqlite_extension))
  return(sqlite_file)
}

import_data = function(myBgeeObject, experimentId = NULL, sampleId = NULL, anatEntityId = NULL, stageId = NULL) {
  user_experiments <- detect_experiments(myBgeeObject, experimentId, sampleId, anatEntityId, stageId)
  missing_experiments <- experiments_to_download(myBgeeObject = myBgeeObject, experimentId = user_experiments, sqlite_file = sqlitePath(myBgeeObject))
  integrate_experiments(myBgeeObject = myBgeeObject, experimentId = missing_experiments, sqlite_file = sqlitePath(myBgeeObject))
} 

# function using annotations to select which experiments the user needs to load all wanted data
detect_experiments = function(myBgeeObject, experimentId = NULL, sampleId = NULL, anatEntityId = NULL, stageId = NULL) {
  annotation <- getAnnotation(myBgeeObject)
  experiments <- annotation$sample.annotation
  # filter list of experiments to download
  if(!is.null(experimentId)) {
    experiments <- experiments[experiments$Experiment.ID %in% experimentId,]
  } 
  if(!is.null(sampleId)) {
    if(myBgeeObject$dataType == "rna_seq") {
      experiments <- experiments[experiments$Library.ID %in% sampleId,]
    } else if(myBgeeObject$dataType == "affymetrix") {
      experiments <- experiments[experiments$Chip.ID %in% sampleId,]
    } else {
      stop("unrecognized datatype : ", myBgeeObject$dataType)
    }
  }
  if(!is.null(anatEntityId)) {
    experiments <- experiments[experiments$Anatomical.entity.ID %in% anatEntityId,]
  }
  if(!is.null(stageId)) {
    experiments <- experiments[experiments$Stage.ID %in% stageId,]
  }
  return(unique(experiments$Experiment.ID))
}

# function that check among experiment wanted by user which one are already stored in the local database.
# return experiments not already present in the database.
experiments_to_download = function(myBgeeObject, experimentId, sqlite_file) {
  if (file.exists(sqlite_file)) {
    conn <- dbConnect(RSQLite::SQLite(), sqlite_file)
    tables <- dbListTables(conn)
    if (isTRUE(myBgeeObject$dataType %in% tables)) {
      exp_queries <- paste0("Select distinct([Experiment.ID]) from ", myBgeeObject$dataType)
      existing_experiments <- dbGetQuery(conn, exp_queries)
      dbDisconnect(conn)
      return(experimentId[(!experimentId %in% existing_experiments$Experiment.ID)])
    }
    dbDisconnect(conn)
  }
  # if local database does not exist return all experiments
  return(experimentId)
}

# function that take care of both downloading and integrating data in the local database
# one local database is created per species (for the moment the package does not contain any
# multispecies functionality)
integrate_experiments = function(myBgeeObject, experimentId, sqlite_file) {
  if (!length(experimentId) == 0) {
    message("downloading data from Bgee FTP...")
    conn <- dbConnect(RSQLite::SQLite(), sqlite_file)
    # If more than 15 experiments have to be downloaded then download all experiments of this 
    # species (in order not to download too many files from the Bgee FTP)
    if (length(experimentId) > 15) {
      message("You tried to download more than 15 experiments, because of that all the Bgee data for this species will be downloaded.")
      myData <- download_and_uncompress_species(myBgeeObject, experimentId)
    } else {
      myData <- download_and_uncompress_experiment(myBgeeObject, experimentId)
    }
    message("Save data in local sqlite database")
    for(i in seq(myData)) {
      dbWriteTable(conn = conn, name = myBgeeObject$dataType, value = myData[i], header=TRUE, 
        append=TRUE, sep="\t")
    }
    unlink(myData)
    dbDisconnect(conn)
  }
}

# function that download all data of one species for the specified datatype from the ftp
# return the name of all uncompressed files
download_and_uncompress_species = function(myBgeeObject, experimentId) {
  message("Downloading all expression data for species ", myBgeeObject$speciesName)
  allExpressionValues <- basename(myBgeeObject$allExperimentsUrl)
  success <- bgee_download_file(url = myBgeeObject$allExperimentsUrl, 
                                destfile = file.path(myBgeeObject$pathToData, allExpressionValues), mode = 'wb')
  if (success != 0){
    stop("ERROR: Download from FTP was not successful.")
  }
  if (grepl(".zip$", allExpressionValues, perl = TRUE)) {
    message("Saved expression data file in", myBgeeObject$pathToData, "folder. Now unzip file...")
    allExperiments <- unzip(file.path(myBgeeObject$pathToData, allExpressionValues), exdir=myBgeeObject$pathToData)
    unlink(file.path(myBgeeObject$pathToData, allExpressionValues))
    tempFiles <- grep(paste0(".*", paste(experimentId, collapse=".*|.*"), ".*"), allExperiments, value = TRUE)
    myData <- unlist(lapply(tempFiles, unzip, exdir=myBgeeObject$pathToData))
    file.remove(allExperiments)
    # at this point remaining .zip files correspond to archive downloaded but not
    # asked by the user. They have to to be removed
    unlink(file.path(myBgeeObject$pathToData, "*.zip"))

  } else if (grepl(".tar.gz$", allExpressionValues, perl = TRUE)) {
    message("Saved expression data file in", myBgeeObject$pathToData, "folder. Now untar file...")
    # When using untar it is only possible to untar OR list files/dir present in a tarball. 
    # It is not possible to do both actions in one line of code
    tempFilesAndDir <- untar(file.path(myBgeeObject$pathToData, allExpressionValues), exdir=myBgeeObject$pathToData)
    tempFilesAndDir <- untar(file.path(myBgeeObject$pathToData, allExpressionValues), exdir=myBgeeObject$pathToData, list = TRUE)
    # keep only path to tarball to integrate to the local database. Other tarballs will not be  integrated
    tempFilesAndDir <- grep(paste0(".*",paste(experimentId,collapse=".*|.*"),".*"), tempFilesAndDir, value = TRUE)
    tempFilesAndDir <- file.path(myBgeeObject$pathToData, tempFilesAndDir)
    tempFiles <- NULL
    for(i in seq(tempFilesAndDir)) {
      # uncompress files and not directories (there was a bug when tar tried to uncompress GTeX experiment)
      # added this verification to be sure there will not be error with other experiments/species
      if (isTRUE(file_test("-f", tempFilesAndDir[i]))) {
        # uncompress expression files
        myData <-untar(tempFilesAndDir[i], exdir=myBgeeObject$pathToData)
        tempFiles <- c(tempFiles, tempFilesAndDir[i])
      } 
    }
    # list all expression files
    myData <- file.path(myBgeeObject$pathToData, unlist(lapply(tempFiles, untar,
      exdir=myBgeeObject$pathToData, list = TRUE)))
    # delete intermediary archives
    unlink(file.path(myBgeeObject$pathToData, allExpressionValues))
    unlink(tempFiles)
    message("Finished uncompress tar files")
    # at this point remaining tar.gz files correspond to tarball downloaded but not
    # asked by the user. They have to to be removed
    unlink(file.path(myBgeeObject$pathToData, "*.tar.gz"))
  }
  return(myData)
}

# function that download a list of experiments from the ftp and return the name of downloaded
# uncompressed files
download_and_uncompress_experiment = function(myBgeeObject, experimentId) {
  message("download experiments")
  tempFiles <- NULL
  for (i in seq(experimentId)) {
    finalExperimentUrl <- gsub("EXPIDPATTERN", experimentId[i], myBgeeObject$experimentUrl)
    tempFile <- file.path(myBgeeObject$pathToData, basename(finalExperimentUrl))
    success <- bgee_download_file(url = finalExperimentUrl, destfile = tempFile, mode = 'wb')
    if (success != 0){
      stop("ERROR: Download from FTP was not successful. Check the experiments present in Bgee with the getAnnotation() function.")
    }
    tempFiles <- c(tempFiles, tempFile)
  }
  message("uncompress processed expression data...")
  if (grepl(".tar.gz$", tempFiles[1], perl = TRUE)) {
    myData <- file.path(myBgeeObject$pathToData, lapply(tempFiles, untar, exdir=myBgeeObject$pathToData))
    myData <- file.path(myBgeeObject$pathToData, unlist(lapply(tempFiles, untar, exdir=myBgeeObject$pathToData, list = TRUE)))
  } else if (grepl(".zip$", tempFiles[1], perl = TRUE)) {
    myData <- unlist(lapply(tempFiles, unzip, exdir=myBgeeObject$pathToData))
  }
  unlink(tempFiles)
  return(myData)
}

# function that query the data in the local database
query_data = function(myBgeeObject, experimentId = NULL, sampleId = NULL, anatEntityId = NULL, stageId = NULL, sqlite_file = sqlitePath(myBgeeObject)) {
  conn <- dbConnect(RSQLite::SQLite(), sqlite_file)
  # generate query
  query <- paste0("SELECT * from ", myBgeeObject$dataType)
  if(! (is.null(experimentId) & is.null(sampleId) & is.null(anatEntityId) & is.null(stageId)) ) {
    query <- paste0(query, " WHERE ")
    if (!is.null(experimentId)) {
      if(length(experimentId) == 1) {
        query <- paste0(query, "[Experiment.ID] = \"",as.character(experimentId), "\"")
      } else {
        query <- paste0(query, "[Experiment.ID] IN (\"",paste(as.character(experimentId), collapse="\", \""), "\")")
      }
    }
    if (!is.null(sampleId)) {
      if (!is.null(experimentId)) {
        query <- paste0(query, " AND ")
      }
      if(myBgeeObject$dataType == "rna_seq") {
        if (length(sampleId) == 1 ) {
          query <- paste0(query, "[Library.ID] = \"", as.character(sampleId), "\"")
        } else {
          query <- paste0(query, "[Library.ID] IN (\"",paste(as.character(sampleId), collapse="\", \""), "\")")
        }
      } else if (myBgeeObject$dataType == "affymetrix") {
        if (length(sampleId) == 1 ) {
          query <- paste0(query, "[Chip.ID] = \"", as.character(sampleId), "\"")
        } else {
          query <- paste0(query, "[Chip.ID] IN (\"",paste(as.character(sampleId), collapse="\", \""), "\")")
        }
      }
    }
    if (!is.null(anatEntityId)) {
      if (!(is.null(experimentId) & is.null(sampleId))) {
        query <- paste0(query, " AND ")
      }
      if(length(anatEntityId) == 1) {
        query <- paste0(query, "[Anatomical.entity.ID] = \"", as.character(anatEntityId), "\"")
      } else {
        query <- paste0(query, "[Anatomical.entity.ID] IN (\"",paste(as.character(anatEntityId), collapse="\", \""), "\")")
      }
    }
    if (!is.null(stageId)) {
      if (!(is.null(experimentId) & is.null(sampleId) & is.null(anatEntityId))) {
        query <- paste0(query, " AND ")
      }
      if(length(stageId) == 1) {
        query <- paste0(query, "[Stage.ID] = \"", as.character(stageId), "\"")
      } else {
        query <- paste0(query, "[Stage.ID] IN (\"",paste(as.character(stageId), collapse="\", \""), "\")")
      }
    }
  }
  message("Load queried data. The query is : ", query)
  result <- dbGetQuery(conn, query)
  dbDisconnect(conn)
  return(result)
}