#' @title Retrieve Bgee RNA-seq or Affymetrix data.
#'
#' @description This function loads the quantitative expression data and presence calls for samples available from Bgee (rna_seq, affymetrix).
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species and data type.
#'
#' @param experimentId An ArrayExpress or GEO accession, e.g., GSE30617. Default is NULL: takes all available experiments for targeted species and data type.
#'
#' @return If experimentId is not specified, returns a list of data frames with data from all experiments for targeted species and data type. If experimentId is specified, returns a data frame with data from this experiment.
#'
#' @author Andrea Komljenovic and Julien Roux.
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'   dataMouse <- getData(bgee)
#'   dataMouseGSE30617 <- getData(bgee, experimentId = "GSE30617")
#' }
#'
#' @export

getData = function(myBgeeObject, experimentId = NULL){

  ## check that fields of Bgee object are not empty
  if (length(myBgeeObject$quantitativeData) == 0 | length(myBgeeObject$experimentUrl) == 0 | length(myBgeeObject$allExperimentsUrl) == 0 | length(myBgeeObject$dataType) == 0 | length(myBgeeObject$pathToData) == 0){
    stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object is valid.")
  }
  ## TO DO: more explicit error when quantitativeData is FALSE.
  ## For example check if length(dataType) > 1 and warn that data download possible for only one data type?
  ## if (length(dataType) == 1 & (dataType == "rna_seq" | dataType == "affymetrix"))
  ## -> maybe specified data type not available for this species
  ## separate testing quantitativeData parameter from others? If quantitativeData == FALSE, annotationUrl will be NULL for example
  ## Do this alos in getAnnotation

  ## Check that download of data is possible for targeted species and data type
  if (myBgeeObject$quantitativeData != TRUE){
    stop("ERROR: downloading the quantitative data is not possible for the species and data type of the input Bgee class object.")
  }

  ## TO DO: code from Bgee.R to integrate here?
  ## check data type availability for chosen species. We only need to verify quantitativeData field (maybe dataType too)
  # if (dataType == "rna_seq" &
  #     allSpecies$RNA_SEQ[allSpecies$ID == speciesId] == FALSE) {
  #   stop(
  #     "ERROR: The specified species name is not among the list of species with RNA-seq data in Bgee release ",
  #     release,
  #     ". See listBgeeSpecies() for details on data types availability for each species."
  #   )
  # } else if (dataType == "affymetrix" &
  #            allSpecies$AFFYMETRIX[allSpecies$ID == speciesId] == FALSE) {
  #   stop(
  #     "ERROR: The specified species name is not among the list of species with Affymetrix microarray data in Bgee release ",
  #     release,
  #     ". See listBgeeSpecies() for details on data types availability for each species."
  #   )
  # }
  #
  ## TO DO: do this also in getAnnotation

  if (length(experimentId) == 0){
    cat(paste0("The experiment is not defined. Hence taking all ", myBgeeObject$dataType, " experiments available for ", myBgeeObject$speciesName, ".\n"))

    ## Get name of data file from URL
    all_expression_values <- basename(myBgeeObject$allExperimentsUrl)

    ## check if RDS file already in cache. If so, skip download step
    if (file.exists(paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_all_experiments_expression_data.rds"))){
      cat(paste0("\nNOTE: expression data file in .rds format was found in the download directory ", myBgeeObject$pathToData,
                 ". Data will not be redownloaded.\n"))
      data_all <- readRDS(file = paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_all_experiments_expression_data.rds"))
    } else {
      cat("Downloading expression data...\n")
      success <- download.file(myBgeeObject$allExperimentsUrl,
                               destfile=file.path(myBgeeObject$pathToData, all_expression_values),
                               mode='wb')
      if (success != 0){
        stop("ERROR: Download from FTP was not successful.")
      }
      cat("Saved expression data file in", myBgeeObject$pathToData, "folder.\n")
      cat("Unzipping file...\n")
      temp.files <- unzip(file.path(myBgeeObject$pathToData, all_expression_values), exdir=myBgeeObject$pathToData)
      mydata <- lapply(temp.files, unzip, exdir=myBgeeObject$pathToData)
      data_all <- lapply(unlist(mydata, rec = TRUE), function(x) as.data.frame(suppressWarnings(fread(x))))

      ## remove spaces in headers
      for (i in 1:length(data_all)){
        names(data_all[[i]]) <- make.names(names(data_all[[i]]))
      }

      cat("Saving all data in .rds file...\n")
      saveRDS(data_all, file = paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_all_experiments_expression_data.rds"))

      ## clean up downloaded files
      file.remove(temp.files)
      file.remove(unlist(mydata, rec = TRUE))
      file.remove(file.path(myBgeeObject$pathToData, all_expression_values))
    }
  } else if( length(experimentId) == 1){
    if (!grepl("^GSE\\d+$|^E-\\w+-\\d+.*$", experimentId, perl = TRUE)){
      stop("The experimentId field needs to be a valid GEO or ArrayExpress accession, e.g., 'GSE30617' or 'E-MEXP-2011'")
    } else {
      ## Since experiment Id is defined, we can now substitute it in the URL
      finalExperimentUrl <- gsub("EXPIDPATTERN", experimentId, myBgeeObject$experimentUrl)
      temp.file <- file.path(myBgeeObject$pathToData, basename(finalExperimentUrl))

      ## check if RDS file already in cache. If so, skip download step
      if (file.exists(paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_", experimentId, "_expression_data.rds"))){
        cat(paste0("\nNOTE: expression data file in .rds format was found in the download directory ", myBgeeObject$pathToData,
                   " for", experimentId, ". Data will not be redownloaded.\n"))
        data_all <- readRDS(paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_", experimentId, "_expression_data.rds"))
      } else {
        cat("Downloading expression data for the experiment", experimentId, "...\n")
        success <- download.file(finalExperimentUrl,
                                 destfile=temp.file,
                                 mode='wb')
        if (success != 0){
          stop("ERROR: Download from FTP was not successful. Check the experiments present in Bgee with the getAnnotation() function.")
        }
        cat("Saved expression data file in", myBgeeObject$pathToData, "folder.\n")
        cat(paste0("Unzipping ", temp.file," file...\n"))
        # Unzipping this file can give one expression data file or multiple ones (if multiple chip types used in experiment)
        mydata <- unzip(temp.file, exdir=myBgeeObject$pathToData)
        data_all <- lapply(mydata, function(x) as.data.frame(fread(x)))

        ## remove spaces in headers
        for (i in 1:length(data_all)){
          names(data_all[[i]]) <- make.names(names(data_all[[i]]))
        }

        if (length(data_all) == 1){
          data_all <- as.data.frame(data_all[[1]])
        }
        cat("Saving all data in .rds file...\n")
        saveRDS(data_all, file = paste0(myBgeeObject$pathToData, "/", myBgeeObject$dataType, "_", experimentId, "_expression_data.rds"))

        ## cleaning up downloaded files
        file.remove(temp.file)
        file.remove(mydata)
      }
    }
  } else {
    stop("Please provide only one experiment accession. If you want to get all data for this species and data type, leave experimentId empty")
  }

  return(data_all)
  cat("Done.")
}
