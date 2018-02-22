#' @title Retrieve Bgee RNA-seq or Affymetrix data.
#'
#' @description This function loads the quantitative expression data and presence calls for samples available from Bgee (rna_seq, affymetrix).
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species and data type.
#'
#' @param experimentId An ArrayExpress or GEO accession, e.g., GSE43721. Default is NULL: takes all available experiments for targeted species and data type.
#'
#' @return If experimentId is not specified, returns a list of data frames with data from all experiments for targeted species and data type. If experimentId is specified, returns a data frame with data from this experiment.
#'
#' @author Andrea Komljenovic and Julien Roux.
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'   dataMouse <- getData(bgee)
#'   dataMouseGSE43721 <- getData(bgee, experimentId = "GSE43721")
#' }
#'
#' @export

getData = function(myBgeeObject, experimentId = NULL){

  ## check that the Bgee object is valid
  if (length(myBgeeObject$quantitativeData) == 0 ){
    stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object was correctly built.")
  } else {
    ## Check that download of quantitative data is possible for targeted species and data type
    ## Try to output error message that helps a bit the user
    if (myBgeeObject$quantitativeData != TRUE){
      if (length(myBgeeObject$dataType) > 1){
        stop("ERROR: downloading quantitative data is only possible if a single data type (\"rna_seq\" or \"affymetrix\") is specified in the input Bgee class object.")
      } else if (length(myBgeeObject$dataType) == 1 & (myBgeeObject$dataType == "rna_seq" | myBgeeObject$dataType == "affymetrix")){
        stop("ERROR: downloading quantitative data is not possible for the species and data type specified in the input Bgee class object. Maybe the specified data type is not available for the targeted species in Bgee? See listBgeeSpecies() for details on data types availability for each species.")
      } else {
        stop("ERROR: downloading quantitative data is not possible for the species and data type specified in the input Bgee class object.")
      }
    } else if (length(myBgeeObject$experimentUrl) == 0 | length(myBgeeObject$allExperimentsUrl) == 0 | length(myBgeeObject$dataType) == 0 | length(myBgeeObject$pathToData) == 0){
      stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object was correctly built.")
    }
  }

  if (length(experimentId) == 0){
    cat(paste0("\nThe experiment is not defined. Hence taking all ", myBgeeObject$dataType, " experiments available for ", myBgeeObject$speciesName, ".\n"))

    ## Get name of data file from URL
    allExpressionValues <- basename(myBgeeObject$allExperimentsUrl)

    ## check if RDS file already in cache. If so, skip download step
    if (file.exists(paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_all_experiments_expression_data.rds"))){
      cat(paste0("\nNOTE: expression data file in .rds format was found in the download directory ", myBgeeObject$pathToData,
                 ". Data will not be redownloaded.\n"))
      allData <- readRDS(file = paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_all_experiments_expression_data.rds"))
    } else {
      cat("\nDownloading expression data...\n")
      success <- download.file(myBgeeObject$allExperimentsUrl,
                               destfile=file.path(myBgeeObject$pathToData, allExpressionValues),
                               mode='wb')
      if (success != 0){
        stop("ERROR: Download from FTP was not successful.")
      }
      if(regexpr("*.zip$", allExpressionValues)){
        cat("\nSaved expression data file in", myBgeeObject$pathToData, "folder. Now unzip file...\n")
        tempFiles <- unzip(file.path(myBgeeObject$pathToData, allExpressionValues), exdir=myBgeeObject$pathToData)
        myData <- lapply(tempFiles, unzip, exdir=myBgeeObject$pathToData)
      }else if(regexpr("*.tar.gz$", allExpressionValues)){
        cat("\nSaved expression data file in", myBgeeObject$pathToData, "folder. Now untar file...\n")
        #using untar it is possible to untar OR list files present in a tarball. It is not possible to do both actions in one line of code
        tempFiles <- untar(file.path(myBgeeObject$pathToData, allExpressionValues), exdir=myBgeeObject$pathToData)
        tempFiles <- untar(file.path(myBgeeObject$pathToData, allExpressionValues), exdir=myBgeeObject$pathToData, list = TRUE)
        myData <- lapply(file.path(myBgeeObject$pathToData, tempFiles), untar, exdir=myBgeeObject$pathToData)
        myData <- lapply(file.path(myBgeeObject$pathToData, tempFiles), untar, exdir=myBgeeObject$pathToData, list = TRUE)
      }else{
        stop("\nThe compressed file can not be unzip nor untar\n")
      }
      allData <- lapply(unlist(myData, recursive = TRUE), function(x) as.data.frame(suppressWarnings(fread(x))))

      ## remove spaces in headers
      for (i in 1:length(allData)){
        names(allData[[i]]) <- make.names(names(allData[[i]]))
      }

      cat("\nSaving all data in .rds file...\n")
      saveRDS(allData, file = paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_all_experiments_expression_data.rds"))

      ## clean up downloaded files
      file.remove(tempFiles)
      file.remove(unlist(myData, recursive = TRUE))
      file.remove(file.path(myBgeeObject$pathToData, allExpressionValues))
    }
  } else if( length(experimentId) == 1){
    if (!grepl("^GSE\\d+$|^E-\\w+-\\d+.*$|^SRP\\d+$", experimentId, perl = TRUE)){
      stop("The experimentId field needs to be a valid GEO, ArrayExpress or SRA accession, e.g., 'GSE43721', 'E-MEXP-2011' or 'SRP044781'")
    } else {
      ## Since experiment Id is defined, we can now substitute it in the URL
      finalExperimentUrl <- gsub("EXPIDPATTERN", experimentId, myBgeeObject$experimentUrl)
      tempFile <- file.path(myBgeeObject$pathToData, basename(finalExperimentUrl))

      ## check if RDS file already in cache. If so, skip download step
      if (file.exists(paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_", experimentId, "_expression_data.rds"))){
        cat(paste0("\nNOTE: expression data file in .rds format was found in the download directory ", myBgeeObject$pathToData,
                   " for ", experimentId, ". Data will not be redownloaded.\n"))
        allData <- readRDS(paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_", experimentId, "_expression_data.rds"))
      } else {
        cat("\nDownloading expression data for the experiment", experimentId, "...\n")
        success <- download.file(finalExperimentUrl,
                                 destfile=tempFile,
                                 mode='wb')
        if (success != 0){
          stop("ERROR: Download from FTP was not successful. Check the experiments present in Bgee with the getAnnotation() function.")
        }
        # Unzipping this file can give one expression data file or multiple ones (if multiple chip types used in experiment)
        if(regexpr("*.zip$",tempFile)){
          cat(paste0("\nSaved expression data file in ", myBgeeObject$pathToData, " folder. Now unzip ", tempFile," file...\n"))
          myData <- unzip(tempFile, exdir=myBgeeObject$pathToData)
        }else if(regexpr("*.tar.gz$",tempFile)){
          cat(paste0("\nSaved expression data file in ", myBgeeObject$pathToData, " folder. Now untar ", tempFile," file...\n"))
          #using untar it is possible to untar OR list files present in a tarball. It is not possible to do both actions in one line of code
          untar(tempFile, exdir=myBgeeObject$pathToData)
          myData <- untar(tempFile, exdir=myBgeeObject$pathToData, list = TRUE)
        }else{
          stop("\nThe file can not be decompressed because it is not a zip nor a tar.gz file\n")
        }
        
        allData <- lapply(myData, function(x) as.data.frame(fread(x)))

        ## remove spaces in headers
        for (i in 1:length(allData)){
          names(allData[[i]]) <- make.names(names(allData[[i]]))
        }

        if (length(allData) == 1){
          allData <- as.data.frame(allData[[1]])
        }
        cat("\nSaving all data in .rds file...\n")
        saveRDS(allData, file = paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_", experimentId, "_expression_data.rds"))

        ## cleaning up downloaded files
        file.remove(tempFile)
        file.remove(myData)
      }
    }
  } else {
    stop("Please provide only one experiment accession. If you want to get all data for this species and data type, leave experimentId empty")
  }

  return(allData)
  cat("\nDone.")
}

