#' @title Retrieve expression of orthologous genes from different species but one single experiment using BgeeDB R package.
#'
#' @description This function allows to load bulk RNA-Seq processed expression values from one experiment with libraries coming
#' from different species. Orthologous genes have the same ``Gene.Family`` identifier
#' 
#' @param bgeeRelease (default current)
#'
#' @param downloadPath path to the directory where expression data will be downloaded
#' 
#' @param experimentId the ID of the experiment for which bulk RNA-Seq expression data have to be retrieved.
#' Can not be NULL. The experiment ID needs to exist and to have libraries in Bgee for provided ``speciesIds``
#' 
#' @param speciesIds species IDs for which we want to retrieve bulk RNA-Seq expression data. Can not be NULL.
#' The ``speciesIds`` needs to exist and to have expression data in Bgee for the specified ``experimentId``.
#' If 1-to-many orthologs have to be retrieved then the reference species has to be provided as ``referenceSpecies`` and should
#' not be added to the list of ``speciesIds``
#' 
#' @param anatEntityIds List of anatomical entity IDs for which expression data have to be retrieved. If NULL then all
#' anatomical entities will be retrieved.
#' 
#' @param removeDownloadFiles Does downloaded files have to be deleted once orthologous genes have been retrieved.
#' Do not delete files by default.
#' 
#' @param onlyOneToOneOrthologs boolean defining whether only one-to-one orthologs have to be retrieved (TRUE) or not (default : FALSE)
#' 
#' @param onlyOneToManyOrthologs boolean defining whether only one-to-many orthologs have to be retrieved (TRUE) or not (default : FALSE)
#' 
#' @param referenceSpecies Correspond to one species ID. Has to be provided only if onlyOneToOne is FALSE. It will then retrieve all ortholog
#' genes having more than one copy in mandatory and optional species.
#' 
#' @param createFile By default the function return the data.frame with expression for all orthologous genes. If ``createFile`` is TRUE then
#' expression data are saved in a file and no data is returned from this function
#' 
#' @param orthologs a list of orthologs to use to retrieve expression. If NULL the list of orthologs is generated using
#' the union of species present in `speciesIds` and `referenceSpeciesId`. It also reuse values of `onlyOneToOneOrthologs`
#' and `onlyOneToManyOrthologs`
#' 
#' @author Julien Wollbrett
#' 
#' @import BgeeDB
#'
#' @examples{
#'   # retrieve expression for 1-to-1 orthologs from 3 species
#'   expression <- retrieveOrthologsExpression(experimentId = "GSE30352", speciesIds = c(9606,10090, 13616), onlyOneToOneOrthologs = TRUE)
#'   # write expression for 1-to-1 orthologs from 3 species in a file
#'   retrieveOrthologsExpression(experimentId = "GSE30352", speciesIds = c(9606,10090, 13616), onlyOneToOneOrthologs = TRUE, createFile = TRUE)
#'   # retrieve expression for genes of 10090 and all of its orthologous genes in 9606 and 13616
#'   expression <- retrieveOrthologsExpression(experimentId = "GSE30352", referenceSpecies = 10090, speciesIds = c(9606, 13616))
#'   # retrieve expression for genes of 10090 and all its duplicated orthologous genes in 9606 and 13616
#'   expression <- retrieveOrthologsExpression(experimentId = "GSE30352", referenceSpecies = 10090, speciesIds = c(9606, 13616), onlyOneToManyOrthologs = TRUE)
#' }
#' 
retrieveOrthologsExpression <- function(bgeeRelease = "current", downloadPath = getwd(),
                                        experimentId = NULL, referenceSpeciesId = NULL, speciesIds = NULL, anatEntityIds = NULL, removeDownloadFiles = FALSE,
                                        onlyOneToOneOrthologs = FALSE, onlyOneToManyOrthologs = FALSE, createFile = FALSE, orthologs = NULL) {
  # we consider that BgeeDB package is already installed... lazy

  # first retrieve ortholog genes for the specified species if not provided
  if (is.null(orthologs)) {
    orthologs <- listOrthologs(mandatorySpecies = speciesIds, referenceSpecies = referenceSpeciesId, onlyOneToOne = onlyOneToOneOrthologs,
                               onlyOneToMany = onlyOneToManyOrthologs, bgeeRelease = bgeeRelease)
  }
  
  # ugly quick and dirty implementation
  gene_families <- unique(orthologs[1])
  
  # now retrieve expression data for all libraries of one species
  all_expression <- NULL
  allSpeciesIds <- c(referenceSpeciesId, speciesIds)
  for (speciesId in allSpeciesIds) {
    bgee_object <- Bgee$new(species = as.character(speciesId), dataType = "rna_seq")
    library_metadata <- getAnnotation(bgee_object)[[1]]
    #only keep metadata from provided experiment
    library_metadata <- library_metadata[library_metadata$Experiment.ID == experimentId,]
    # now filter on anatomical entities
    if (!is.null(anatEntityIds)) {
      library_metadata <- library_metadata[library_metadata$Anatomical.entity.ID %in% anatEntityIds,]
    }
    if (nrow(library_metadata) == 0) {
      stop("No library available for species ", speciesId, ". Please remove that species or update the experimentId, ",
           "anatEntityIds or bgeeRelease")
    }
    expression_data <- getData(sampleId = unique(library_metadata$Library.ID), myBgeeObject = bgee_object)
    expression_data$Species.ID <- speciesId
    orthologs_expression_data <- expression_data[expression_data$Gene.ID %in% orthologs[,as.character(speciesId)],]
    #implementation with lapply (slow so commented it)
    # which_f <- function(x, orthologs, speciesId, gene_families) {return(which(gene_families == orthologs[which(orthologs[,as.character(speciesId)] == x),1]))}
    # orthologs_expression_data$Gene.Family <- lapply(X = orthologs_expression_data$Gene.ID, FUN = which_f, orthologs = orthologs, speciesId = speciesId, gene_families = gene_families)
    # naive implementation with for loop... faster than using lapply
    orthologs_expression_data$Gene.Family <- NA
    for (row_number in seq(nrow(orthologs_expression_data))) {
      geneId <- orthologs_expression_data$Gene.ID[row_number]
      orthologs_row_indice <- min(which(orthologs[,as.character(speciesId)] == geneId))
      gene_family_row_indice <- which(gene_families == orthologs[orthologs_row_indice,1])
      orthologs_expression_data$Gene.Family[row_number] <- gene_family_row_indice
    }
    all_expression <- rbind(all_expression, orthologs_expression_data)
  }

  if (createFile) {
    output_file <- file.path(downloadPath, paste0("expression_", paste(allSpeciesIds, collapse = "_"), ".tsv"))
    message("writing expression data in the file ", output_file)
    write.table(x = all_expression, file = output_file,
                sep = "\t", quote = F, row.names = F, col.names = T)
  } else {
    message("properly retrieved Bgee expression data")
    return(all_expression)
  }
}

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
#' @param referenceSpecies Correspond to one species ID. Has to be provided only if onlyOneToOne is FALSE.
#' It will then retrieve all ortholog genes having more than one copy in mandatory and optional species.
#' 
#' @param onlyOneToOne boolean defining whether only one-to-one orthologs have to be retrieved (TRUE) or not (default : FALSE)
#' 
#' @param onlyOneToMany boolean defining whether only one-to-many orthologs have to be retrieved (TRUE) or not (default : FALSE)
#' 
#' @param removeDownloadedFiles Does downloaded files have to be deleted once orthologous genes have been retrieved.
#' Do not delete files by default (FALSE).
#' 
#' @param createFile By default the function return the data.frame of all orthologous genes. If ``createFile`` is TRUE then the
#' orthologous genes are saved in a file and no data is returned from this function
#' 
#' @author Julien Wollbrett
#'
#' @examples{
#'   # 1-to-1 orthologs
#'   orthologs <- listOrthologs(mandatorySpecies = c(9606,10090,13616), onlyOneToOne = TRUE)
#'   # 1-to-1 orthologs saved in a file
#'   listOrthologs(mandatorySpecies = c(9606,10090,7955,7227,13616), onlyOneToOne = TRUE, createFile = TRUE)
#'   # all orthologs (1-to-1 and 1-to-many) from 10090, 9606 and 13616 with potential duplication in 9606 or 13616
#'   orthologs <- listOrthologs(referenceSpecies = 10090, mandatorySpecies = c(9606,13616))
#'   # 1-to-many orthologs from 10090, 9606 and 13616 with duplication in 9606 and 13616
#'   orthologs <- listOrthologs(referenceSpecies = 10090, mandatorySpecies = c(9606,13616), onlyOneToMany = TRUE)
#' }
#' 
listOrthologs <- function(bgeeRelease = "current", downloadPath = getwd(), 
                          mandatorySpecies = NULL, optionalSpecies = NULL, referenceSpecies = NULL, onlyOneToOne = FALSE,
                          onlyOneToMany = FALSE, removeDownloadedFiles = FALSE, createFile = FALSE) {
  
  # TODO check all arguments combination and throw error if not an expected one (e.g both onlyOneToOne and onlyOneToMany)
  if (!is.null(referenceSpecies) && onlyOneToOne) {
    stop("please provide a reference species only when 1-to-many orthologs have to be retrieved")
  }
  if (onlyOneToMany && onlyOneToOne) {
    stop("please select only one argument out of onlyOneToOne and onlyOneToMany")
  }
  if (is.null(referenceSpecies) && !onlyOneToOne) {
    stop("please provide a reference species if 1-to-many orthologs have to be retrieved")
  }
  
  # first draft of the implementation lots of future functionalities to implement....
  if (!is.null(optionalSpecies)) {stop("optionalSpecies not yet implemented")}
  
  ftp_url <- "https://www.bgee.org/ftp/RELEASE_VERSION/homologous_genes/OMA_orthologs.zip"
  orthologs_dir <- "bgeeOrthologs"
  archive_file <- basename(ftp_url)
  download_dir <- file.path(downloadPath, orthologs_dir)
  unzipped_files <- NULL
  if (!dir.exists(download_dir)) {
    message("dir does not exist")
    dir.create(download_dir)
    downloaded_archive <- file.path(download_dir, archive_file)
    ftp_url <- gsub(pattern = "RELEASE_VERSION", replacement = bgeeRelease, x = ftp_url)
    # for now we only provide one zip archive containing all orthologs files
    file_path <- download.file(url = ftp_url, destfile = downloaded_archive)
    unzipped_files <- unzip(zipfile = downloaded_archive, exdir = download_dir)
    unzipped_files <- list.files(path = download_dir, include.dirs = FALSE, full.names = TRUE, recursive = TRUE,
                                 pattern = "orthologs_.*.csv")
  } else {
    message("Bgee orthologs already downloaded")
    unzipped_files <- list.files(path = download_dir, include.dirs = FALSE, full.names = TRUE, recursive = TRUE,
                                 pattern = "orthologs_.*.csv")
  }
  
  # now files are downloaded and orthologs retrieval can starts
  #for now we consider retrieval of all orthologs from mandatory species
  retrieved_orthologs_files <- NULL
  retrieved_species <- NULL
  all_orthologs <- NULL
  for (file_path in unzipped_files) {
    speciesIds <- unlist(regmatches(x = basename(file_path), m = gregexpr(pattern = "[0-9]+", basename(file_path))))
    if (!is.null(referenceSpecies) && referenceSpecies %in% as.numeric(speciesIds) && 
        (as.numeric(speciesIds[1]) %in% mandatorySpecies || as.numeric(speciesIds[2]) %in% mandatorySpecies)) {
      orthologs <- read.table(file = file_path, header = TRUE, sep = ",", quote = "")
      if (speciesIds[1] == as.character(referenceSpecies)) {
        orthologs <- orthologs[,c(1,2)]
        colnames(orthologs) <- c(referenceSpecies, speciesIds[2])
      } else {
        orthologs <- orthologs[,c(2,1)]
        colnames(orthologs) <- c(referenceSpecies, speciesIds[1])
      }
      duplicated_species2 <- unique(orthologs[duplicated(orthologs[2]),][2])
      orthologs <- orthologs[! orthologs[,2] %in% duplicated_species2[[1]],]
      if (onlyOneToMany) {
        duplicated_species1 <- unique(orthologs[duplicated(orthologs[1]),][1])
        orthologs <- orthologs[ orthologs[,1] %in% duplicated_species1[[1]],]
      }
      if (as.character(referenceSpecies) %in% retrieved_species) {
        all_orthologs <- merge(x = all_orthologs, y = orthologs, by = as.character(referenceSpecies))
      } else {
        all_orthologs <- orthologs
      }
      retrieved_species <- unique(c(retrieved_species, speciesIds))
    } else if (is.null(referenceSpecies) && as.numeric(speciesIds[1]) %in% mandatorySpecies &&
               as.numeric(speciesIds[2] %in% mandatorySpecies)) {
      retrieved_orthologs_files <- c(retrieved_orthologs_files, file_path)
      orthologs <- read.table(file = file_path, header = TRUE, sep = ",", quote = "")[,1:2]
      # only one-to-one requested so we remove all duplicated genes
      duplicated_species2 <- unique(orthologs[duplicated(orthologs[2]),][2])
      duplicated_species1 <- unique(orthologs[duplicated(orthologs[1]),][1])
      orthologs <- orthologs[! orthologs[,1] %in% duplicated_species1[[1]],]
      orthologs <- orthologs[! orthologs[,2] %in% duplicated_species2[[1]],]
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
  if (removeDownloadedFiles) {
    file.remove(download_dir, recursive = TRUE)
  }
  if (createFile) {
    output_file <- file.path(downloadPath, paste0("orthologs_", paste(retrieved_species, collapse = "_"), ".tsv"))
    message("writing orthologous genes in the file ", output_file)
    write.table(x = all_orthologs, file = output_file,
                sep = "\t", quote = F, row.names = F, col.names = T)
  } else {
    message("properly retrieved Bgee orthologs")
    return(all_orthologs)
  }
}

#' @title Format data.frame for 1-to-N orthologs
#'
#' @description This function take as input the output of the listOrthologs() function and format it to have one line
#' per gene family (i.e one line per gene from the reference species having orthologs in all other species). The resulting
#' dataframe has 1 + (referenceSpeciesId * n) columns. column names correspond to the the speciesId followed by a suffix. For
#' instance for 1-to-2 orthologs where the reference species is the gar and the other species is D. rerio, then columns will be
#' 7918 (gar species ID), 7955_1 (first ortholog in D. rerio), 7955_2 (2nd ortholog in D. rerio) 
#'
#' @param orthologs the output of the listOrthologs() function
#'
#' @param n the exact number of times the gene is duplicated in geneDuplicatedSpecies species. Correspond to the n
#' in 1-to-n orthologs
#' 
#' @param referenceSpecies the species for which orthologs have to be retrieved. correspond to the 1 in 1-to-n orthologs
#'
#' @param geneDuplicatedSpecies Correspond to a list of speciesId for which duplicated orthologs have been retrieved
#' 
#' @author Julien Wollbrett
#'
#' @examples{
#'   # 1-to-many orthologs in Gar and zebrafish with gar being the reference species
#    orthologs_one_to_many <- listOrthologs(referenceSpecies = 7918, mandatorySpecies = c(7955), onlyOneToMany  = TRUE)
#'   # reformat and only keep 1-to-2 orthologs
#'   one_to_two_ortholgos <- reformat_orthologs_one_to_n(orthologs = orthologs_one_to_many, n = 2, referenceSpecies = 7918,
#'   geneDuplicateSpecies = 7955)
#' }
#' 

formatOrthologsOneToN <- function(orthologs = NULL, n = 2, referenceSpecies = NULL,
                                        geneDuplicatedSpecies = NULL) {
  one_to_n_orhtologs <- as.data.frame(matrix(nrow = 0, ncol = 1+(n * length(geneDuplicatedSpecies))))
  header <- referenceSpecies
  for(speciesId in geneDuplicatedSpecies) {
    for (i in seq_len(n)) {
      header <- c(header, paste0(speciesId, "_", i))
    }
  }
  colnames(one_to_n_orhtologs) <- header
  for(reference_geneId in unique(orthologs[,c(as.character(referenceSpecies))])) {
    current_one_to_n_geneIds <- reference_geneId
    for(speciesId in geneDuplicatedSpecies) {
      gene_n_species <- orthologs[orthologs[,c(as.character(referenceSpecies))] == reference_geneId,]
      gene_n_species <- unique(gene_n_species[,c(as.character(referenceSpecies), as.character(speciesId))])
      if(nrow(gene_n_species) == n) {
        merged_gene_n_species <- NULL
        for (i in seq_len(n)) {
          current_one_to_n_geneIds <- c(current_one_to_n_geneIds, gene_n_species[i,2])
        }
      }
    }
    current_one_to_n_geneIds <- as.data.frame(t(current_one_to_n_geneIds))
    rownames(current_one_to_n_geneIds) <- NULL
    if(ncol(current_one_to_n_geneIds) == length(header)) {
      colnames(current_one_to_n_geneIds) <- header
      one_to_n_orhtologs <- rbind(one_to_n_orhtologs, current_one_to_n_geneIds)
    }
  }
  return(one_to_n_orhtologs)
}