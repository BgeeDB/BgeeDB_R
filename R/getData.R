########################
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

getData = function(..., experimentId = NULL){

  if (length(experiment.id) == 0){
    cat(paste0("The experiment is not defined. Hence taking all ", dataType, " experiments available for ", speciesName, ".\n"))

    all_expression_values <- basename(allexperimentsUrl)

    ## check if RDS file already in cache. If so, skip download step
    if (file.exists(paste0(pathToData, "/", dataType, "_all_experiments_expression_data.rds"))){
      cat(paste0("NOTE: expression data file in .rds format was found in the download directory ", pathToData,
                 ". Data will not be redownloaded.\n"))
      data_all <- readRDS(file = paste0(pathToData, "/", dataType, "_all_experiments_expression_data.rds"))
    } else {
      cat("Downloading expression data...\n")
      success <- download.file(allexperimentsUrl,
                               destfile=file.path(pathToData, all_expression_values),
                               mode='wb')
      if (success != 0){
        stop("ERROR: Download from FTP was not successful.")
      }
      cat("Saved expression data file in", pathToData, "folder.\n")
      cat("Unzipping file...\n")
      temp.files <- unzip(file.path(pathToData, all_expression_values), exdir=pathToData)
      mydata <- lapply(temp.files, unzip, exdir=pathToData)
      data_all <- lapply(unlist(mydata, rec = TRUE), function(x) as.data.frame(suppressWarnings(fread(x))))

      ## remove spaces in headers
      for (i in 1:length(data_all)){
        names(data_all[[i]]) <- make.names(names(data_all[[i]]))
      }

      cat("Saving all data in .rds file...\n")
      saveRDS(data_all, file = paste0(pathToData, "/", dataType, "_all_experiments_expression_data.rds"))

      ## cleaning up downloaded files
      file.remove(temp.files)
      file.remove(unlist(mydata, rec = TRUE))
      file.remove(file.path(pathToData, all_expression_values))
    }
  } else if( length(experiment.id) == 1){
    if (!grepl("^GSE\\d+$|^E-\\w+-\\d+.*$", experiment.id, perl = TRUE)){
      stop("The experiment.id field needs to be a GEO or ArrayExpress accession, e.g., 'GSE30617' or 'E-MEXP-2011'")
    } else {
      cat("Downloading expression data for the experiment", experiment.id, "\n")

      ## Since experiment Id is defined now, we can substitute it in the URL
      experimentUrl <<- gsub("EXPIDPATTERN", experiment.id, experimentUrl)
      temp.file <- file.path(pathToData, basename(experimentUrl))

      ## check if RDS file already in cache. If so, skip download step
      if (file.exists(paste0(pathToData, "/", dataType, "_", experiment.id, "_expression_data.rds"))){
        cat(paste0("NOTE: expression data file in .rds format was found in the download directory ", pathToData,
                   " for", experiment.id, ". Data will not be redownloaded.\n"))
        data_all <- readRDS(paste0(pathToData, "/", dataType, "_", experiment.id, "_expression_data.rds"))
      } else {
        cat("Downloading expression data...\n")
        success <- download.file(experimentUrl,
                                 destfile=temp.file,
                                 mode='wb')
        if (success != 0){
          stop("ERROR: Download from FTP was not successful. Check the experiments present in Bgee with the get_annotation() function.")
        }
        cat("Saved expression data file in", pathToData, "folder.\n")
        cat(paste0("Unzipping ", temp.file," file...\n"))
        # Unzipping this file can give one expression data file or multiple ones (if multiple chip types used in experiment)
        mydata <- unzip(temp.file, exdir=pathToData)
        data_all <- lapply(mydata, function(x) as.data.frame(fread(x)))

        ## remove spaces in headers
        for (i in 1:length(data_all)){
          names(data_all[[i]]) <- make.names(names(data_all[[i]]))
        }

        if (length(data_all) == 1){
          data_all <- as.data.frame(data_all[[1]])
        }
        cat("Saving all data in .rds file...\n")
        saveRDS(data_all, file = paste0(pathToData, "/", dataType, "_", experiment.id, "_expression_data.rds"))

        ## cleaning up downloaded files
        file.remove(temp.file)
        file.remove(mydata)
      }
    }
  } else {
    stop("Please provide only one experiment ID. If you want to get all data for this species and data type, leave experiment.id empty")
  }

  return(data_all)
  cat("Done.")
}
