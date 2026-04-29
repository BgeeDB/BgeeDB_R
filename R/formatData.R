#' @title Format RNA-seq or Affymetrix data downloaded from Bgee.
#'
#' @description This function formats the data downloaded with the getSampleRawData() function into an object of the Bioconductor "expressionSet" Class.
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species and data type.
#'
#' @param data A list of data frames including data from multiple experiments, or a data frame including data from a single experiment.
#'
#' @param stats A character indicating what expression values should be used in the formatted data expressionSet object matrix.
#'  \itemize{
#'    \item{"rpkm" for bulk RNA-seq (Bgee release 13.2 and before)}
#'    \item{"fpkm" for bulk RNA-seq (Bgee release 14.0 to 14.2)}
#'    \item{"counts" for RNA-seq}
#'    \item{"tpm" for bulk RNA-seq or full length single cell RNA-seq (Bgee release 14 and above)}
#'    \item{"cpm" for droplet based single cell RNA-seq (Bgee release 15.2 and above)}
#'    \item{"intensities" for Affymetrix microarrays}
#'  }
#'
#' @param callType A character indicating whether intensities should be displayed only for present (i.e., expressed) genes, present high quality genes, or all genes (default).
#'  \itemize{
#'    \item{"present"}
#'    \item{"present high quality"}
#'    \item{"all"}
#'  }
#'
#' @return If data was a list of data frames from multiple experiments, returns a list of ExpressionSet objects. If data was a data frame from a single experiment, returns an ExpressionSet object.
#'
#' @author Andrea Komljenovic and Julien Roux.
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'   dataMouseGSE43721 <- getData(bgee, experimentId = "GSE43721")
#'   dataMouseGSE43721.tpm <- formatData(bgee,
#'                                        dataMouseGSE43721,
#'                                        callType = "present",
#'                                        stats = "tpm")
#' }
#'
#' @importFrom dplyr %>%
#' @importFrom tidyr spread_
#' @export

formatData = function(myBgeeObject, data, stats = NULL, callType = "all"){

  ## check that fields of Bgee object are not empty
  if (length(myBgeeObject$quantitativeData) == 0 | length(myBgeeObject$dataType) == 0){
    stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object is valid.")
  }

  ## Check that download of data is possible for targeted species and data type
  if (myBgeeObject$quantitativeData != TRUE){
    stop("ERROR: formatting quantitative data is not possible for the species and data type of the input Bgee class object.")
  }

  if (length(stats) == 0){
    if (myBgeeObject$dataType == "affymetrix"){
      cat("\nWARNING: stats parameter set to \"intensities\" for the next steps.\n")
      stats <- "intensities"
    } else if (myBgeeObject$dataType %in% c("rna_seq", "sc_full_length")) {
      stop("Please specify the stats parameters. Should be set to \"counts\", \"tpm\" (For Bgee 14.0 and above), \"rpkm\" (for Bgee 13.2) or \"fpkm\" (between Bgee 14.0 and 14.2)")
    } else if (myBgeeObject$dataType %in% c("sc_droplet_based")) {
      stop("Please specify the stats parameters. Should be set to \"counts\" or \"cpm\"")
    }
  # manage stats error for all datatypes and Bgee releases
  } else if (myBgeeObject$dataType == "affymetrix" & stats != "intensities"){
    cat("\nWARNING: For Affymetrix microarray data, stats parameter can only be set to \"intensities\". This will be used for the next steps.\n")
    stats <- "intensities"
  } else if (myBgeeObject$dataType == "rna_seq" & compareVersion(gsub("_", ".", myBgeeObject$release), "13.2") <= 0 & !(stats %in% c('rpkm', 'counts'))){
    stop("Choose whether data formatting should create a matrix of RPKMs or read counts, with stats option set as \"rpkm\" or \"counts\"")
  } else if (myBgeeObject$dataType == "rna_seq" & compareVersion(gsub("_", ".", myBgeeObject$release), "14.0") >= 0 
      & compareVersion(gsub("_", ".", myBgeeObject$release), "14.2") <= 0 & !(stats %in% c('fpkm', 'counts', 'tpm'))){
    stop("Choose whether data formatting should create a matrix of FPKMs, TPMs or read counts, with stats option set as \"fpkm\", \"tpm\" or \"counts\"")
  } else if (myBgeeObject$dataType == "rna_seq" & compareVersion(gsub("_", ".", myBgeeObject$release), "14.2") > 0 & !(stats %in% c('counts', 'tpm'))){
    stop("Choose whether data formatting should create a matrix of TPMs or read counts, with stats option set as \"tpm\" or \"counts\"")
  } else if (myBgeeObject$dataType == "sc_full_length" & !(stats %in% c('counts', 'tpm'))){
    stop("Choose whether data formatting should create a matrix of TPMs or read counts, with stats option set as \"tpm\" or \"counts\"")
  } else if (myBgeeObject$dataType == "sc_droplet_based" & !(stats %in% c('counts', 'cpm'))){
    stop("Choose whether data formatting should create a matrix of CPMs or read counts, with stats option set as \"cpm\" or \"counts\"")
  }

  if(!(callType %in% c('present','present high quality','all'))){
    stop("Choose between displaying intensities for present genes, present high quality genes or all genes, e.g., 'present', 'present high quality', 'all' ")
  }

  if(length(data) == 1) data[[1]] else data

  cat("\nExtracting expression data matrix...\n")
  if(stats  == "rpkm"){
    columns <- c("Library.ID", "Gene.ID", "RPKM")
    expr <- .extract.data(data, columns, callType)
  } else if(stats  == "fpkm"){
    columns <- c("Library.ID", "Gene.ID", "FPKM")
    expr <- .extract.data(data, columns, callType)
  } else if(stats  == "tpm"){
    columns <- c("Library.ID", "Gene.ID", "TPM")
    expr <- .extract.data(data, columns, callType)
  } else if (stats == "counts"){
    columns <- c("Library.ID", "Gene.ID", "Read.count")
    expr <- .extract.data(data, columns, callType)
  } else if(stats  == "cpm"){
    columns <- c("Library.ID", "Gene.ID", "CPM")
    expr <- .extract.data(data, columns, callType)
  } else {
    columns <- c("Chip.ID", "Probeset.ID", "Log.of.normalized.signal.intensity", "Gene.ID")
    expr <- .extract.data(data, columns, callType)
  }
  ## one data matrix
  if(is.data.frame(expr$assayData)){

    ## sort objects to have samples in the same order
    expr$assayData <- expr$assayData[, order(names(expr$assayData))]
    expr$pheno <- expr$pheno[order(sampleNames(expr$pheno))]

    ## create ExpressionSet object
    eset <- new('ExpressionSet',
                exprs=as.matrix(expr$assayData),
                phenoData=expr$pheno,
                featureData=expr$features)
  } else if(is.list(expr$assayData)){
    ## multiple data matrices
    eset <- mapply(function(x,y,z){
      ## sort objects to have samples in the same order
      x <- x[, order(names(x))]
      y <- y[order(sampleNames(y))]

      ## create ExpressionSet object
      new('ExpressionSet',
          exprs=as.matrix(x),
          phenoData=y,
          featureData=z)
    },
    expr$assayData, expr$pheno, expr$features)
  }
  return(eset)
}

## Take in unformatted data downloaded from Bgee and outputs a list of expression matrices, phenotype annotations, feature annotations, and calls.
.extract.data <- function(data, columns, callType){
  ## if multiple data frames (multiple experiments or chips)
  if(class(data) == "list"){
    calls <- lapply(data, function(x) .calling(x, callType, columns[3]))
    expr <- lapply(calls,
                   function(x) {
                     ## subset the data to keep relevant columns
                     xt <- x[, columns[1:3]]
                     ## from sample and feature columns, create a matrix with features as rows and samples as columns
                     xtt <- xt %>% spread_(columns[1], columns[3])
                     rownames(xtt) <- xtt[, columns[2]]
                     ## Remove feature column to keep only data
                     xtt[,-1, drop = FALSE]
                   }
    )
    cat("\nExtracting features information...\n")
    features <- mapply(.extract.data.feature, calls, expr, rep(list(columns), times=length(calls)))
    cat("\nExtracting samples information...\n")
    phenos <- mapply(.extract.data.pheno, calls, rep(list(columns[1]), times=length(calls)))
  } else {
    ## if only a single dataframe
    calls <- .calling(data, callType, columns[3])
    ## subset the data to keep relevant columns
    xt <- calls[, columns[1:3]]
    xtt <- xt %>% spread_(columns[1], columns[3])
    rownames(xtt) <- xtt[, columns[2]]
    ## Remove feature column to keep only data
    expr <- xtt[,-1, drop = FALSE]

    cat("\nExtracting features information...\n")
    features <- .extract.data.feature( calls, expr, columns)
    cat("\nExtracting samples information...\n")
    phenos <- .extract.data.pheno( calls, columns[1])
  }
  return(list(assayData = expr, pheno = phenos, features = features, calls = calls))
}

# Extract feature data (probesets or genes)
.extract.data.feature <- function(calls, expr, columns){
  ## RNA-seq
  if(length(columns) == 3){
    fdata <- calls[match(rownames(expr), calls[, columns[2]]), columns[2], drop = FALSE]
  }
  ## Affymetrix, 4 columns
  else if(length(columns) == 4){
    fdata <- calls[match(rownames(expr), calls[, columns[2]]), columns[c(2,4)], drop = FALSE]
  }
  rownames(fdata) <- fdata[, columns[2]]
  fdata <- as(fdata, "AnnotatedDataFrame")
  return(fdata)
}

# Extract annotation of samples
.extract.data.pheno <- function(calls, column){
  phdata <- calls[, c(column, "Anatomical.entity.ID", "Anatomical.entity.name", "Stage.ID", "Stage.name")]
  phdata <- phdata[!duplicated(phdata[, column]), ]
  rownames(phdata) <- phdata[, column]
  phdata <- as.data.frame(phdata)
  metadata <- data.frame(labelDescription=colnames(phdata),
                         row.names=colnames(phdata))
  phenodata <- new("AnnotatedDataFrame",
                   data=phdata,
                   varMetadata=metadata
  )
  return(phenodata)
}

.calling <- function(x, callType, column){
  ## check data type
  if(callType == "present"){
    cat("  Keeping only present genes.\n")
    x[(x$Detection.flag == "absent"), column] <- NA
  } else if (callType == "present high quality"){
    cat("  Keeping only present high quality genes.\n")
    x[which(x$Detection.flag == "absent" | x$Detection.quality == "poor quality"), column] <- NA
  }
  return(x)
}
