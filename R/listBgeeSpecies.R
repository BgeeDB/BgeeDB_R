

.listDirectories <- function(url, version = FALSE) {
	tmpcon <- textConnection(getURL(url), "r")
	tx <- read.table(tmpcon)
	close(tmpcon)
	ifelse(version, tx1 <- paste0(levels(tx$V6), tx$V8)[1], tx1 <- as.character(tx[,ncol(tx)]))
	return(tx1)
}


#' @title List of species in Bgee database
#' @description Takes a URL and returns a list of available genomes for different platforms in Bgee.
#' @param ... an empty parameter
#'
#'
#' @return A list of available RNA-seq and Affymetrix datasets
#' @import data.table RCurl
#' @export

listBgeeSpecies <- function(...){
	url_rnaseq <-  "ftp://ftp.bgee.org/current/download/processed_expr_values/rna_seq/"
	url_affymetrix <-  "ftp://ftp.bgee.org/current/download/processed_expr_values/affymetrix/"
	cat("Last update:" , .listDirectories(url_rnaseq, version = TRUE), "\n")
	return(list(rna_seq = .listDirectories(url_rnaseq), affymetrix = .listDirectories(url_affymetrix)))

}
