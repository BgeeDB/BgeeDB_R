#' @title Retrieve Bgee RNA-seq or Affymetrix data.
#'
#' @description This function loads the quantitative expression data and presence calls 
#' for samples available from Bgee (rna_seq, affymetrix, sc_full_length).
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
#' @param cellTypeId Filter specific to single cell datatype (sc_full_length) allowing to specify 
#' one or more cell type IDs from the UBERON ontology (http://uberon.github.io/). Default is 
#' NULL: takes all available cell types for targeted species and data type. Available for Bgee 15.0 and after
#' 
#' @param sex Filter allowing to specify one or more sexes. Default is 
#' NULL: takes all available sexes for targeted species and data type. Available for Bgee 15.0 and after
#' 
#' @param strain Filter allowing to specify one or more strains. Default is 
#' NULL: takes all available strains for targeted species and data type. Available for Bgee 15.0 and after
#' 
#' @param withDescendantAnatEntities Allows to filter on the selected anatEntityId and all its descendants.
#' This functionality is available for Bgee 15.0 release and after
#'
#' @param withDescendantStages Allows to filter on the selected stageId and all its descendants.
#' This functionality is available for Bgee 15.0 release and after
#'
#' @param withDescendantCellTypes Allows to filter on the selected cellTypeId and all its descendants.
#' This functionality is available for Bgee 15.0 release and after
#'
#' @return Return a dataframe containing all Bgee processed expression data from the selected species 
#' and datatype using specified filters with operator AND.
#'
#' @author Julien Wollbrett, Andrea Komljenovic and Julien Roux.
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'   dataMouseGSE43721 <- getData(bgee, experimentId = "GSE43721")
#'   dataMouseVariousFilters <- getData(bgee, experimentId = c("GSE43721", "GSE36026"), 
#'                              anatEntityId = c("UBERON:0002107", "UBERON:0000956", "UBERON:0002048"))
#' }
#'
#' @import RSQLite
#' @importFrom R.utils gunzip
#' @export
#' 
getData <- function(myBgeeObject, experimentId = NULL, sampleId = NULL, 
    anatEntityId = NULL, stageId = NULL, cellTypeId = NULL, sex = NULL, strain = NULL, 
    withDescendantAnatEntities = FALSE, withDescendantStages = FALSE, 
    withDescendantCellTypes = FALSE) {
  check_object(myBgeeObject)
  # write a warning message if user tried to retrieve ontology terms with descendants for
  # a bgee release where this functionality was not yet implemented
  compare_version_number <- gsub("_", ".", myBgeeObject$release)
  if ((withDescendantAnatEntities | withDescendantStages | withDescendantCellTypes) &
      compareVersion(a = compare_version_number , b = "15.0") < 0) {
    message("withDescendant functionality is available only for Bgee 15.0",
            " release and after. Will not retrieved descendant of selected parameters.")
  }
  if (withDescendantAnatEntities & compareVersion(a = compare_version_number , b = "15.0") >= 0) {
    if(is.null(anatEntityId)) {
      warning("No anatomical entity was provided. Not possible to filter on descendant anatomical entities.")
    }
    anatEntityId <- c(anatEntityId, getDescendantAnatEntities(bgee = myBgeeObject, ids = anatEntityId))
  }
  if (withDescendantStages & compareVersion(a = compare_version_number , b = "15.0") >= 0) {
    if(is.null(stageId)) {
      warning("No developmental stage was provided. Not possible to filter on descendant developmental stages.")
    }
    stageId <- c(stageId, getDescendantStages(bgee = myBgeeObject, ids = stageId))
  }
  if (withDescendantCellTypes & compareVersion(a = compare_version_number , b = "15.0") >= 0) {
    if(is.null(cellTypeId)) {
      warning("No cell type was provided. Not possible to filter on descendant cell types.")
    }
    cellTypeId <- c(cellTypeId, getDescendantCellTypes(bgee = myBgeeObject, ids = cellTypeId))
  }
  check_condition_parameters(myBgeeObject = myBgeeObject, anatEntityId = anatEntityId, 
    stageId = stageId, cellTypeId = cellTypeId, sex = sex, strain = strain)
  import_data(myBgeeObject = myBgeeObject, experimentId = experimentId, sampleId = sampleId, 
    anatEntityId = anatEntityId, stageId = stageId, cellTypeId = cellTypeId, sex = sex, 
    strain = strain)
  return(query_data(myBgeeObject = myBgeeObject, experimentId = experimentId, sampleId = sampleId, 
    anatEntityId = anatEntityId, stageId = stageId, cellTypeId = cellTypeId, sex = sex, 
    strain = strain))
}

check_object = function(myBgeeObject, experimentId = NULL){
  ## check that the Bgee object is valid
  if (length(myBgeeObject$quantitativeData) == 0 ){
    stop("ERROR: there seems to be a problem with the input Bgee class object, some fields ", 
         "are empty. Please check that the object was correctly built.")
  } else {
    ## Check that download of quantitative data is possible for targeted species and data type
    ## Try to output error message that helps a bit the user
    if (myBgeeObject$quantitativeData != TRUE){
      if (length(myBgeeObject$dataType) > 1){
        stop("ERROR: downloading quantitative data is only possible if a single data ", 
             "type (\"rna_seq\" or \"affymetrix\") is specified in the input Bgee class object.")
      } else if (length(myBgeeObject$dataType) == 1 & !(myBgeeObject$dataType %in% c('rna_seq','affymetrix'))){
        stop("ERROR: downloading quantitative data is not possible for the species and data ", 
        "type specified in the input Bgee class object. Maybe the specified data type is not available ", 
        "for the targeted species in Bgee? See listBgeeSpecies() for details on data types availability for each species.")
      } else {
        stop("ERROR: downloading quantitative data is not possible for the species and data type ", 
             "specified in the input Bgee class object.")
      }
    } else if (length(myBgeeObject$experimentUrl) == 0 | length(myBgeeObject$allExperimentsUrl) == 0 
        | length(myBgeeObject$dataType) == 0 | length(myBgeeObject$pathToData) == 0){
      stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are ", 
      "empty. Please check that the object was correctly built.")
    }
  }
}

check_condition_parameters = function(myBgeeObject, anatEntityId, stageId, 
    cellTypeId, sex, strain){
  ## check that the condition parameters queried are compatible with Bgee release
  ## selected
  if(compareVersion(a = gsub("_", ".", myBgeeObject$release), b = "15.0") < 0) {
    if (!is.null(cellTypeId) | !is.null(sex) | !is.null(strain)) {
      stop("ERROR: cellTypeId, sex, and strain can be filtered only for Bgee 15.0",
      " release and after.")
    }
  }
}

## extract the global path to the sqliteDB
sqlitePath <- function(myBgeeObject){
  sqlite_file <- file.path(myBgeeObject$pathToData, paste0(myBgeeObject$speciesName, 
    myBgeeObject$sqlite_extension))
  return(sqlite_file)
}

import_data = function(myBgeeObject, experimentId = NULL, sampleId = NULL, 
    anatEntityId = NULL, stageId = NULL, cellTypeId = NULL, sex = NULL, strain = NULL) {
  user_experiments <- detect_experiments(myBgeeObject, experimentId, sampleId, anatEntityId, 
    stageId, cellTypeId, sex, strain)
  # stop code if no experiments correspond to user criteria
  if(length(user_experiments) == 0) {
    stop("No data correspond to selected parameters. Please look at the annotation to ",
         "select available data.")
  }
  missing_experiments <- experiments_to_download(myBgeeObject = myBgeeObject, 
    experimentId = user_experiments, sqlite_file = sqlitePath(myBgeeObject))
  integrate_experiments(myBgeeObject = myBgeeObject, experimentId = missing_experiments, 
    sqlite_file = sqlitePath(myBgeeObject))
} 

# function using annotations to select which experiments the user needs to load all wanted data
detect_experiments = function(myBgeeObject, experimentId = NULL, sampleId = NULL, 
    anatEntityId = NULL, stageId = NULL, cellTypeId = NULL, sex = NULL, strain = NULL) {
  annotation <- suppressMessages(getAnnotation(myBgeeObject))
  experiments <- annotation$sample.annotation
  # filter list of experiments to download
  if(!is.null(experimentId)) {
    experiments <- experiments[experiments$Experiment.ID %in% experimentId,]
  } 
  if(!is.null(sampleId)) {
    if(myBgeeObject$dataType == "rna_seq" | myBgeeObject$dataType == "sc_full_length") {
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
  if(!is.null(sex)) {
    # annotation files contain NA values. dplyr:na_if(sex, "NA") function allows to change
    # queried "NA" values to NA values
    experiments <- experiments[experiments$Sex %in% na_if(sex, "NA"),]
  }
  if(!is.null(strain)) {
    experiments <- experiments[experiments$Strain %in% strain,]
  }
  if(!is.null(cellTypeId)) {
    if (myBgeeObject$dataType == "sc_full_length") {
      experiments <- experiments[experiments$Cell.type.ID %in% cellTypeId,]
    } else {
      stop("Can only filter on cell type ID when single cell datatype (sc_full_length) is selected.")
    }
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
    # If more than 15 experiments have to be downloaded then download all experiments of this 
    # species (in order not to download too many files from the Bgee FTP)
    if (length(experimentId) > 15) {
      message("You tried to download more than 15 experiments, because of that all the ", 
        "Bgee data for this species will be downloaded.")
      myData <- download_and_uncompress_species(myBgeeObject, experimentId)
    } else {
      myData <- download_and_uncompress_experiment(myBgeeObject, experimentId)
    }
    message("Save data in local sqlite database")
    for(i in seq(myData)) {
      conn <- dbConnect(RSQLite::SQLite(), sqlite_file)
      dbWriteTable(conn = conn, name = myBgeeObject$dataType, 
        value = myData[i], header=TRUE, append=TRUE, sep="\t")
      dbDisconnect(conn)
      # those lines are mandatory as Rank and pValue can have NA values that
      # are transformed to 0 by dbQuery() function.
    }
    if(compareVersion(a = gsub("_", ".", myBgeeObject$release), b = "15.0") >= 0) {
      conn <- dbConnect(RSQLite::SQLite(), sqlite_file)
      updatePvalues <- dbExecute(conn, paste0("UPDATE ",myBgeeObject$dataType," set [pValue] = NULL 
        where [pValue] = \"NA\""))
      dbDisconnect(conn)
      if(myBgeeObject$dataType == "rna_seq" | myBgeeObject$dataType == "sc_full_length") {
        conn <- dbConnect(RSQLite::SQLite(), sqlite_file)
        updateRanks <- dbExecute(conn, paste0("UPDATE ", myBgeeObject$dataType, " set [Rank] = NULL 
          where [Rank] = \"NA\""))
        dbDisconnect(conn)
      }
    }
    unlink(myData)
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
  # compressed files for Bgee 13.2
  if (grepl(".zip$", allExpressionValues, perl = TRUE)) {
    message("Saved expression data file in", myBgeeObject$pathToData, 
      "folder. Now unzip file...")
    allExperiments <- unzip(file.path(myBgeeObject$pathToData, allExpressionValues), 
      exdir=myBgeeObject$pathToData)
    unlink(file.path(myBgeeObject$pathToData, allExpressionValues))
    tempFiles <- grep(paste0(".*", paste(experimentId, collapse=".*|.*"), ".*"), 
      allExperiments, value = TRUE)
    myData <- unlist(lapply(tempFiles, unzip, exdir=myBgeeObject$pathToData))
    file.remove(allExperiments)
  # compressed tarball for bgee14.0 and after
  } else if (grepl(".tar.gz$", allExpressionValues, perl = TRUE)) {
    message("Saved expression data file in", myBgeeObject$pathToData, "folder. Now untar file...")
    # When using untar it is only possible to untar OR list files/dir present in a tarball. 
    # It is not possible to do both actions in one line of code
    tempFilesAndDir <- untar(file.path(myBgeeObject$pathToData, allExpressionValues), 
      exdir=myBgeeObject$pathToData)
    tempFilesAndDir <- untar(file.path(myBgeeObject$pathToData, allExpressionValues), 
      exdir=myBgeeObject$pathToData, list = TRUE)
    # keep only path to tarball to integrate to the local database. Other tarballs will 
    # not be  integrated
    tempFilesAndDir <- grep(paste0(".*",paste(experimentId,collapse=".*|.*"),".*"), 
      tempFilesAndDir, value = TRUE)
    tempFilesAndDir <- file.path(myBgeeObject$pathToData, tempFilesAndDir)
    tempFiles <- NULL
    myData <- NULL
    for(i in seq(tempFilesAndDir)) {
      # uncompress files and not directories (there was a bug when tar tried to uncompress 
      # GTeX experiment)
      # added this verification to be sure there will not be error with other experiments/species
      if (isTRUE(file_test("-f", tempFilesAndDir[i]))) {
        # uncompress expression files for Bgee 14
        if (grepl(".tar.gz$", tempFilesAndDir[i], perl = TRUE)) {
          currentData <-untar(tempFilesAndDir[i], exdir=myBgeeObject$pathToData)
          # list all expression files
          currentData <- file.path(myBgeeObject$pathToData, untar(tempFilesAndDir[i], 
            exdir=myBgeeObject$pathToData, list = TRUE))
        # uncompress expression files for Bgee 15 and after
        } else if (grepl(".tsv.gz$", tempFilesAndDir[i], perl = TRUE)) {
          fileNameUncompressed <- gsub(".gz", "", basename(tempFilesAndDir[i]))
          currentData <- unlist(gunzip(tempFilesAndDir[i], remove=FALSE, overwrite = TRUE,
            destname = file.path(myBgeeObject$pathToData, fileNameUncompressed)))
        }
        tempFiles <- c(tempFiles, tempFilesAndDir[i])
        myData <- c(myData, currentData)
      } 
    }
    # delete intermediary archives
    unlink(file.path(myBgeeObject$pathToData, allExpressionValues))
    unlink(dirname(tempFiles[1]), recursive = TRUE)
    message("Finished uncompress tar files")
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
      stop("ERROR: Download from FTP was not successful. Check the experiments present ", 
      "in Bgee with the getAnnotation() function.")
    }
    tempFiles <- c(tempFiles, tempFile)
  }
  message("uncompress processed expression data...")
  if (grepl(".tar.gz$", tempFiles[1], perl = TRUE)) {
    myData <- file.path(myBgeeObject$pathToData, lapply(tempFiles, untar, 
      exdir=myBgeeObject$pathToData))
    myData <- file.path(myBgeeObject$pathToData, unlist(lapply(tempFiles, untar, 
      exdir=myBgeeObject$pathToData, list = TRUE)))
  } else if (grepl(".zip$", tempFiles[1], perl = TRUE)) {
    myData <- unlist(lapply(tempFiles, unzip, exdir=myBgeeObject$pathToData))
  } else if (grepl(".tsv.gz$", tempFiles[1], perl = TRUE)) {
    myData <- unlist(lapply(tempFiles, gunzip, remove=FALSE, overwrite = TRUE))
  }
  unlink(tempFiles)
  return(myData)
}

# function that query the data in the local database
query_data = function(myBgeeObject, experimentId = NULL, sampleId = NULL, anatEntityId = NULL, 
    stageId = NULL, cellTypeId = NULL, sex = NULL, strain = NULL, 
    sqlite_file = sqlitePath(myBgeeObject)) {
  conn <- dbConnect(RSQLite::SQLite(), sqlite_file)
  # generate query
  query <- paste0("SELECT * from ", myBgeeObject$dataType)
  if(! (is.null(experimentId) & is.null(sampleId) & is.null(anatEntityId) & is.null(stageId)
        & is.null(sex) & is.null(strain) & is.null(cellTypeId)) ) {
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
        query <- paste0(query, "[Stage.ID] IN (\"",paste(as.character(stageId), 
          collapse="\", \""), "\")")
      }
    }
    # All following filters are available for Bgee 15.0 and after
    if(compareVersion(a = gsub("_", ".", myBgeeObject$release), b = "15.0") >= 0) {
      if (!is.null(cellTypeId)) {
        if (!(is.null(experimentId) & is.null(sampleId) & is.null(anatEntityId) 
            & is.null(stageId))) {
          query <- paste0(query, " AND ")
        }
        if(length(cellTypeId) == 1) {
          query <- paste0(query, "[Cell.type.ID] = \"", as.character(cellTypeId), "\"")
        } else {
          query <- paste0(query, "[Cell.type.ID] IN (\"",paste(as.character(cellTypeId), 
            collapse="\", \""), "\")")
        }
      }
      if (!is.null(sex)) {
        if (!(is.null(experimentId) & is.null(sampleId) & is.null(anatEntityId) 
            & is.null(stageId) & is.null(cellTypeId))) {
          query <- paste0(query, " AND ")
        }
        if(length(sex) == 1) {
          query <- paste0(query, "[Sex] = \"", as.character(sex), "\"")
        } else {
          query <- paste0(query, "[Sex] IN (\"",paste(as.character(sex), 
            collapse="\", \""), "\")")
        }
      }
      if (!is.null(strain)) {
        if (!(is.null(experimentId) & is.null(sampleId) & is.null(anatEntityId) 
              & is.null(stageId) & is.null(cellTypeId) & is.null(sex))) {
          query <- paste0(query, " AND ")
        }
        # rsqlite does not allow to remove double quotes when inserting data in the local database.
        # strain data are surrounded by " in the tsv files. These quotes allows to have
        # strain with a name containing special characters.
        # Following the SQL syntax, it is mandatory to provide 2 double quotes in order
        # escape existing double quotes in the database.
        if(length(strain) == 1) {
          query <- paste0(query, "[Strain] = \"\"\"", as.character(strain), "\"\"\"")
        } else {
          query <- paste0(query, "[Strain] IN (\"\"\"",paste(as.character(strain), 
            collapse="\"\"\", \"\"\""), "\"\"\")")
        }
      }
    }
  }
  message("Load queried data. The query is : ", query)
  result <- dbGetQuery(conn, query)
  dbDisconnect(conn)
  return(result)
}

getDescendantAnatEntities <- function (bgee, ids) {
  return(getDescendant(bgee = bgee, ids = ids, conditionParam = "anatEntities"))
}

getDescendantCellTypes <- function (bgee, ids) {
  return(getDescendant(bgee = bgee, ids = ids, conditionParam = "anatEntities"))
}

getDescendantStages <- function (bgee, ids) {
  return(getDescendant(bgee = bgee, ids = ids, conditionParam = "stages"))
}

getDescendant <- function (bgee, ids, conditionParam) {
  myUrl <- paste0(bgee$topAnatUrl,
  "?page=r_package&action=COND_PARAM&ENTITIES&species_id=SPECIES&",
  "propagation=DESCENDANTS&attr_list=ID&display_type=tsv")
  myUrl <- gsub("SPECIES", bgee$speciesId, myUrl, perl = FALSE)
  if (conditionParam == "anatEntities") {
    myUrl <- gsub("COND_PARAM", "get_propagation_anat_entity", myUrl, perl = TRUE)
    myUrl <- gsub("ENTITIES", paste0("anat_entity_id=",
      paste(ids, collapse = "&anat_entity_id=")), myUrl, perl = TRUE)
  } else if (conditionParam == "stages") {
    myUrl <- gsub("COND_PARAM", "get_propagation_dev_stage", myUrl, perl = TRUE)
    myUrl <- gsub("ENTITIES", paste0("stage_id=",
      paste(ids, collapse = "&stage_id=")), myUrl, perl = TRUE)
  }
  destFile <- file.path(bgee$pathToData, "descendants.tsv")
  success <- bgee_download_file(url = myUrl, destfile = destFile)
  descendants <- read.table(destFile, header = TRUE, sep = "\t")
  # retrieve annotation in order to keep only descendants present in the annotation
  annotation <- suppressMessages(getAnnotation(bgee)$sample.annotation)
  if (conditionParam == "anatEntities") {
    present <- unique(annotation$Anatomical.entity.ID)
  } else if (conditionParam == "stages") {
    present <- unique(annotation$Stage.ID)
  }
  # Do not use column name because it is not the same depending on the
  # queried condition parameter
  wanted_descendants <- descendants[descendants[,1] %in% present,]
  unlink(destFile)
  return(wanted_descendants)
}