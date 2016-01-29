


listDirectories <- function(url, version = FALSE) {
  # Takes a URL and returns a character vector of filenames.
  #
  # Args:
  #   url: URL is a ftp site of processed data from Bgee database for RNA seq or Affymetrix.
  #   version: If TRUE, prints last update of the database; if not, not. Default is FALSE.
  #
  # Returns:
  #   A character vector of filenames

	require(RCurl)
  # message(url)
	tmpcon <- textConnection(getURL(url), "r")
	tx <- read.table(tmpcon)
	close(tmpcon)
	ifelse(version, tx1 <- paste0(levels(tx$V6), tx$V8)[1], tx1 <- as.character(tx[,ncol(tx)]))
	return(tx1)
}



# what is available in bgee currently
listBgeeSpecies <- function(...){
	# Takes a URL and returns a list of available genomes for different platforms in Bgee.
  #
  # Args: No args needed
  #   
 	#
  # Returns:
  #   List of available RNA-seq and Affymetrix datasets
	
	url_rnaseq <-  "ftp://ftp.bgee.org/current/download/processed_expr_values/rna_seq/"
	url_affymetrix <-  "ftp://ftp.bgee.org/current/download/processed_expr_values/affymetrix/"
	cat("Last update:" , listDirectories(url_rnaseq, version = TRUE), "\n")	
	return(list(rna_seq = listDirectories(url_rnaseq), affymetrix = listDirectories(url_affymetrix)))

}