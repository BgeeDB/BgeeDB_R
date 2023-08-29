#' @title Retrieve orthologs from Bgee
#'
#' @description This function allows to load all ortholog genes based on criteria provided
#' as arguments of this function. It is for instance possible to retrieve only one-to-one orthologs
#' or one-to-many orthologs. Please read the description of the arguments to properly tune the orthologs
#' to retrieve 
#'
#' @param bgeeRelease (default current)
#'
#' @param downloadPath path to the directory where orthologs will be downloaded
#' 
#' @param mandatorySpecies species for which it is mandatory to find orthologs.
#' 
#' @param optionalSpecies species for which we want to add orthologs if available. Will still retrieve
#' the orthologous genes from mandatory species even if no ortholog in these optionalSpecies
#'
#' @param referenceSpecies Correspond to one species ID. Has to be provided only if onlyOneToMany is TRUE.
#' It will then retrieve all ortholog genes having more than one copy in mandatory and optional species.
#' 
#' @param onlyOneToOne boolean defining whether only one-to-one orthologs have to be retrieved (TRUE) or not (default : FALSE)
#' 
#' @param onlyOneToMany boolean defining whether only one-to-many orthologs have to be retrieved (TRUE) or not (default : FALSE)
#' 
#' @param removeDownloadFiles Does downloaded files have to be deleted once orthologous genes have been retrieved.
#' Do not delete files by default.
#' 
#' @author Julien Wollbrett
#'
#' @examples{
#'   listOrthologs(mandatorySpecies = c(9606,10090,7955,7227,13616), onlyOneToOne = TRUE)
#' }
#' 
listOrthologs <- function(bgeeRelease = "current", downloadPath = getwd(), 
                    mandatorySpecies = NULL, optionalSpecies = NULL, referenceSpecies = NULL, onlyOneToOne = FALSE,
                    onlyOneToMany = FALSE, removeDownloadFiles = FALSE) {
  
  # TODO check all arguments combination and throw error if not an expected one (e.g both onlyOneToOne and onlyOneToMany)
  
  # first draft of the implementation lots of future functionalities to implement....
  message(optionalSpecies)
  if (!is.null(optionalSpecies)) {stop("optionalSpecies not yet implemented")}
  if (!is.null(referenceSpecies)) {stop("referenceSpecies not yet implemented")}
  if (onlyOneToMany) {stop("onlyOneToMany not yet implemented")}
  if (removeDownloadFiles) {stop("removeDownloadFiles not yet implemented")}
  if (onlyOneToMany) {stop("onlyOneToMany not yet implemented")}
  
  ftp_url <- "https://www.bgee.org/ftp/RELEASE_VERSION/homologous_genes/OMA_orthologs.zip"
  orthologs_dir <- "bgeeOrthologs"
  archive_file <- basename(ftp_url)
  download_dir <- file.path(downloadPath, orthologs_dir)
  if(!dir.exists(download_dir)) {dir.create(download_dir)}
  downloaded_archive <- file.path(downloadPath, orthologs_dir, archive_file)
  
  ftp_url <- gsub(pattern = "RELEASE_VERSION", replacement = bgeeRelease, x = ftp_url)
  
  # for now we only provide one zip archive containing all orthologs files
  file_path <- download.file(url = ftp_url, destfile = downloaded_archive)
  unzipped_files <- unzip(zipfile = downloaded_archive, exdir = download_dir)

  # now files are downloaded and orthologs retrieval can starts
  #for now we consider retrieval of all orthologs from mandatory species
  retrieved_orthologs_files <- NULL
  retrieved_species <- NULL
  all_orthologs <- NULL
  for (file_path in unzipped_files) {
    speciesIds <- unlist(regmatches(x = basename(file_path), m = gregexpr(pattern = "[0-9]+", basename(file_path))))
    if ((as.numeric(speciesIds[1]) %in% mandatorySpecies) & (as.numeric(speciesIds[2] %in% mandatorySpecies))) {
      retrieved_orthologs_files <- c(retrieved_orthologs_files, file_path)
      orthologs <- read.table(file = file_path, header = TRUE, sep = ",", quote = "")[,1:2]
      # if only one-to-one then remove all duplicated genes
      if (onlyOneToOne) {
        duplicated_species2 <- unique(orthologs[duplicated(orthologs[2]),][2])
        duplicated_species1 <- unique(orthologs[duplicated(orthologs[1]),][1])
        orthologs <- orthologs[! orthologs[,1] %in% duplicated_species1[[1]],]
        orthologs <- orthologs[! orthologs[,2] %in% duplicated_species2[[1]],]
      }
      colnames(orthologs) <- c(speciesIds[1], speciesIds[2])
      if ((speciesIds[1] %in% retrieved_species) & (speciesIds[2] %in% retrieved_species)) {
        all_orthologs <- merge(x = all_orthologs, y = orthologs, by = speciesIds)
      } else if (speciesIds[1] %in% retrieved_species) {
        all_orthologs <- merge(x = all_orthologs, y = orthologs, by = speciesIds[1])
      } else if (speciesIds[2] %in% retrieved_species) {
        all_orthologs <- merge(x = all_orthologs, y = orthologs, by = speciesIds[2])
      } else {
        all_orthologs <- orthologs
      }
      retrieved_species <- unique(c(retrieved_species, speciesIds))
    }
  }
  write.table(x = all_orthologs, file = file.path(download_dir, paste0("orthologs_", paste(retrieved_species, collapse = "_"), ".tsv")),
    sep = "\t", quote = F, row.names = F, col.names = T)
}

listOrthologs(mandatorySpecies = c(9606,10090,7955,7227,13616), onlyOneToOne = TRUE)
