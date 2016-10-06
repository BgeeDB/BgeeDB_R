#' @title Retrieve Bgee experiments annotation for targeted species and data type.
#'
#' @description This function loads the annotation of experiments and samples of quantitative expression datasets (rna_seq, affymetrix) that are available from Bgee.
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species and data type.

#' @return A list of two elements, including a data frame of the annotation of experiments for chosen species (field "experiment_annotation") and a data frame of the annotation of chips/libraries from these experiments (field "sample_annotation").
#'
#' @author Andrea Komljenovic and Julien Roux.
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'   myAnnotation <- getAnnotation(bgee)
#' }
#'
#' @export

getAnnotation = function(myBgeeObject){

  ## check that fields of Bgee object are not empty
  if (length(myBgeeObject$quantitativeData) == 0 | length(myBgeeObject$annotationUrl) == 0 | length(myBgeeObject$dataType) == 0 | length(myBgeeObject$pathToData) == 0){
   stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object is valid.")
  }

  ## Check that download of data is possible for targeted species and data type
  if (myBgeeObject$quantitativeData != TRUE){
    stop("ERROR: downloading the annotation files is not possible for the species and data type of the input Bgee class object.")
  }

  ## Get name of annotation file from URL
  annotation.file <- basename(myBgeeObject$annotationUrl)

  ## To get the names of experiment and sample files, we start from the annotation.file
  if (myBgeeObject$dataType == "affymetrix"){
    annotation.experiments <- gsub("_chips", "", annotation.file)
  } else if (myBgeeObject$dataType == "rna_seq"){
    annotation.experiments <- gsub("_libraries", "", annotation.file)
  }
  annotation.experiments <- gsub("zip", "tsv", annotation.experiments)
  annotation.samples     <- gsub("_experiments", "", annotation.file)
  annotation.samples     <- gsub("zip", "tsv", annotation.samples)

  ## Check if file is already in cache. If so, skip download step
  if (file.exists(file.path(myBgeeObject$pathToData, annotation.experiments)) && file.exists(file.path(myBgeeObject$pathToData, annotation.samples))){
    cat(paste0("\nNOTE: annotation files for this species were found in the download directory ", myBgeeObject$pathToData,
               ". Data will not be redownloaded.\n"))
  } else {
    cat("Downloading annotation files...\n")
    success <- download.file(myBgeeObject$annotationUrl,
                             destfile=file.path(myBgeeObject$pathToData, annotation.file),
                             mode='wb')
    if (success != 0){
      stop("ERROR: Download from FTP was not successful.")
    }
    unzip(file.path(myBgeeObject$pathToData, annotation.file), exdir=myBgeeObject$pathToData)
    cat("Saved annotation files in", myBgeeObject$pathToData, "folder.\n")
    ## Clean directory and remove .zip file
    file.remove(file.path(myBgeeObject$pathToData, annotation.file))
    ## Test if extracted files are OK
    if (!(file.exists(file.path(myBgeeObject$pathToData, annotation.experiments)) & file.exists(file.path(myBgeeObject$pathToData, annotation.samples)))){
      stop("ERROR: extraction of annotation files from downloaded zip file went wrong.")
    }
  }

  ## Read the 2 annotation files
  myanno <- list(sample_annotation=as.data.frame(fread(file.path(myBgeeObject$pathToData, annotation.samples))),
                 experiment_annotation=as.data.frame(fread(file.path(myBgeeObject$pathToData, annotation.experiments)))
  )
  ## remove spaces in headers
  for (i in 1:length(myanno)){
    names(myanno[[i]]) <- make.names(names(myanno[[i]]))
  }
  return(myanno)
}
