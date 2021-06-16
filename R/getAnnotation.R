#' @title Retrieve Bgee experiments annotation for targeted species and data type.
#'
#' @description This function loads the annotation of experiments and samples of quantitative expression datasets (rna_seq, affymetrix, sc_full_length) that are available from Bgee.
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species and data type.

#' @return A list of two elements, including a data frame of the annotation of experiments for chosen species (field "experiment.annotation") and a data frame of the annotation of chips/libraries from these experiments (field "sample.annotation").
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

  ## check that the Bgee object is valid
  if (length(myBgeeObject$quantitativeData) == 0 ){
    stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object was correctly built.")
  } else {
    ## Check that download of quantitative data is possible for targeted species and data type
    ## Try to output error message that helps a bit the user
    if (myBgeeObject$quantitativeData != TRUE){
      if (length(myBgeeObject$dataType) > 1){
        stop("ERROR: downloading annotation files is only possible if a single data type (\"rna_seq\", \"sc_full_length\" or \"affymetrix\") is specified in the input Bgee class object.")
      } else if (length(myBgeeObject$dataType) == 1 & (myBgeeObject$dataType == "rna_seq" | myBgeeObject$dataType == "sc_full_length" | myBgeeObject$dataType == "affymetrix")){
        stop("ERROR: downloading annotation files is not possible for the species and data type specified in the input Bgee class object. Maybe the specified data type is not available for the targeted species in Bgee? See listBgeeSpecies() for details on data types availability for each species.")
      } else {
        stop("ERROR: downloading annotation files is not possible for the species and data type specified in the input Bgee class object.")
      }
    } else if (length(myBgeeObject$annotationUrl) == 0 | length(myBgeeObject$dataType) == 0 | length(myBgeeObject$pathToData) == 0){
      stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object was correctly built.")
    }
  }

  ## Get name of annotation file from the URL field of Bgee object
  annotationFile <- basename(myBgeeObject$annotationUrl)
  ## To get the names of experiment and sample files, we start from the annotationFile
  if (myBgeeObject$dataType == "affymetrix"){
    annotationExperiments <- gsub("_chips", "", annotationFile)
  } else if (myBgeeObject$dataType == "rna_seq" | myBgeeObject$dataType == "sc_full_length"){
    annotationExperiments <- gsub("_libraries", "", annotationFile)
  }
  annotationSamples     <- gsub("_experiments", "", annotationFile)
  if(grepl(".zip$", annotationExperiments, perl = TRUE)){
    annotationExperiments <- gsub("zip", "tsv", annotationExperiments)
    annotationSamples     <- gsub("zip", "tsv", annotationSamples)
  }else if(grepl(".tar.gz$", annotationExperiments, perl = TRUE)){
    annotationExperiments <- gsub("tar.gz", "tsv", annotationExperiments)
    annotationSamples     <- gsub("tar.gz", "tsv", annotationSamples)
  }else{
    stop("\nThe annotation file can not be uncompressed because it is not a zip nor a tar.gz file\n")
  }
  ## Check if file is already in cache. If so, skip download step
  if (file.exists(file.path(myBgeeObject$pathToData, annotationExperiments)) && file.exists(file.path(myBgeeObject$pathToData, annotationSamples))){
    cat(paste0("\nNOTE: annotation files for this species were found in the download directory ", myBgeeObject$pathToData,
               ". Data will not be redownloaded.\n"))
  } else {
    success <- bgee_download_file(url = myBgeeObject$annotationUrl,
                             destfile = file.path(myBgeeObject$pathToData, annotationFile), mode='wb')
    if (success != 0){
      stop("ERROR: Download from FTP was not successful.")
    }
    if(grepl(".zip$", annotationFile, perl = TRUE)){
      myData <- unzip(file.path(myBgeeObject$pathToData, annotationFile), exdir=myBgeeObject$pathToData)
    }else if(grepl(".tar.gz$", annotationFile, perl = TRUE)){
      untar(file.path(myBgeeObject$pathToData, annotationFile), exdir=myBgeeObject$pathToData)
      myData <- untar(file.path(myBgeeObject$pathToData, annotationFile), exdir=myBgeeObject$pathToData, list = TRUE)
      myData <- file.path(myBgeeObject$pathToData, myData)
      # move annotation files from "experiments repository" to pathToData repository in order to be consistent with zip uncompression
      file.rename(from = myData,  to = file.path(myBgeeObject$pathToData, basename(myData)))
      # delete the "experiment repository" folder
      unlink(dirname(myData), recursive = TRUE)
    }else{
      stop("\nThe annotation file can not be uncompressed because it is not a zip nor a tar.gz file\n")
    }
    cat("\nSaved annotation files in", myBgeeObject$pathToData, "folder.\n")
    ## Clean directory and remove .zip file
    file.remove(file.path(myBgeeObject$pathToData, annotationFile))
    ## Test if extracted files are OK
    if (!(file.exists(file.path(myBgeeObject$pathToData, annotationExperiments)) & file.exists(file.path(myBgeeObject$pathToData, annotationSamples)))){
      stop("ERROR: extraction of annotation files from downloaded zip file went wrong.")
    }
  }

  ## Read the 2 annotation files
  myAnnotation <- list(sample.annotation=as.data.frame(fread(file.path(myBgeeObject$pathToData, annotationSamples))),
                       experiment.annotation=as.data.frame(fread(file.path(myBgeeObject$pathToData, annotationExperiments)))
  )
  ## remove spaces in headers
  names(myAnnotation$experiment.annotation) <- make.names(names(myAnnotation$experiment.annotation))
  names(myAnnotation$sample.annotation) <- make.names(names(myAnnotation$sample.annotation))
  ## Double dots are now present in some column names, due to column headers in FTP file such as "Min. read length"
  ## and "Max. read length". These are replaced by a single dot.
  names(myAnnotation$sample.annotation) <- gsub("\\.\\.", ".", names(myAnnotation$sample.annotation))

  return(myAnnotation)
}
