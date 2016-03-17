
.listDirectories <- function(url, version = FALSE) {
	tmpcon <- textConnection(getURL(url), "r")
	tx <- read.table(tmpcon)
	close(tmpcon)
	ifelse(version, tx1 <- paste0(levels(tx$V6), tx$V8)[1], tx1 <- as.character(tx[,ncol(tx)]))
	return(tx1)
}


#' @title List species with Affymetrix or RNA-seq data in the Bgee database
#' @description Returns a list of available genomes for different platforms in Bgee.
#' @param ... an empty parameter
#'
#'
#' @return A list of species with available Affymetrix or RNA-seq data in the Bgee
#'
#' @author Andrea Komljenovic \email{andrea.komljenovic at unil.ch}.
#' @import data.table RCurl
#' @export

listBgeeSpecies <- function(...){
	url_rnaseq <-  "ftp://ftp.bgee.org/current/download/processed_expr_values/rna_seq/"
	url_affymetrix <-  "ftp://ftp.bgee.org/current/download/processed_expr_values/affymetrix/"
	cat("Last update:" , .listDirectories(url_rnaseq, version = TRUE), "\n")
	return(list(rna_seq = .listDirectories(url_rnaseq), affymetrix = .listDirectories(url_affymetrix)))

}

################################
## list species
## http://localhost:8080/?page=dao&action=org.bgee.model.dao.api.species.SpeciesDAO.getAllSpecies&display_type=tsv&attr_list=ID&attr_list=GENUS&attr_list=SPECIES_NAME&attr_list=COMMON_NAME
## TO DO: this should be the basis of the listBgeeSpecies function (cleaner, more maintenable, etc)
## TO DO: should output data type too
################################

