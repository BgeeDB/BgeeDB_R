########################
#' @title Formats results of the enrichment test on anatomical structures.
#'
#' @description This function loads the results from the topGO test and creates an output table with organ names,
#' fold enrichment and FDR. Data are sorted by p-value and only terms below the specified FDR cutoff are included.
#'
#' @param topAnatData A list produced by the function loadTopAnatData().
#'
#' @param topAnatObject An object produced by the function topAnat().
#'
#' @param results A result object, produced by the runtest() function of topGO.
#'
#' @param cutoff An FDR cutoff between 0 and 1. Only terms with FDR lower than this cutoff are included.
#' Default is 1, meaning that all terms are included.
#'
#' @return A data frame with significantly enriched anatomical structures, sorted by p-value.
#'
#' @author Julien Roux
#'
#' @examples{
#'  ## Launch topGO test on data loaded from Bgee
#'   myTopAnatData <- loadTopAnatData(species = "Mus_musculus", datatype = "rna_seq")
#'   geneList <- as.factor(c(rep(0, times=90), rep(1, times=10)))
#'   names(geneList) <- c("ENSMUSG00000064370", "ENSMUSG00000064368", "ENSMUSG00000064367",
#'                     "ENSMUSG00000064363", "ENSMUSG00000065947", "ENSMUSG00000064360",
#'                     "ENSMUSG00000064358", "ENSMUSG00000064357", "ENSMUSG00000064356",
#'                     "ENSMUSG00000064354", "ENSMUSG00000064351", "ENSMUSG00000064345",
#'                     "ENSMUSG00000064341", "ENSMUSG00000029757", "ENSMUSG00000079941",
#'                     "ENSMUSG00000053367", "ENSMUSG00000016626", "ENSMUSG00000037816",
#'                     "ENSMUSG00000036781", "ENSMUSG00000022519", "ENSMUSG00000079606",
#'                     "ENSMUSG00000068966", "ENSMUSG00000038608", "ENSMUSG00000047473",
#'                     "ENSMUSG00000038542", "ENSMUSG00000025386", "ENSMUSG00000028145",
#'                     "ENSMUSG00000024816", "ENSMUSG00000020978", "ENSMUSG00000055373",
#'                     "ENSMUSG00000038155", "ENSMUSG00000046408", "ENSMUSG00000030032",
#'                     "ENSMUSG00000042249", "ENSMUSG00000071909", "ENSMUSG00000039670",
#'                     "ENSMUSG00000032501", "ENSMUSG00000054252", "ENSMUSG00000068071",
#'                     "ENSMUSG00000067578", "ENSMUSG00000074892", "ENSMUSG00000027905",
#'                     "ENSMUSG00000058216", "ENSMUSG00000078754", "ENSMUSG00000062101",
#'                     "ENSMUSG00000043633", "ENSMUSG00000071350", "ENSMUSG00000021639",
#'                     "ENSMUSG00000059113", "ENSMUSG00000049115", "ENSMUSG00000053310",
#'                     "ENSMUSG00000043832", "ENSMUSG00000063767", "ENSMUSG00000026775",
#'                     "ENSMUSG00000038537", "ENSMUSG00000078716", "ENSMUSG00000096820",
#'                     "ENSMUSG00000075089", "ENSMUSG00000049971", "ENSMUSG00000014303",
#'                     "ENSMUSG00000056054", "ENSMUSG00000033082", "ENSMUSG00000020801",
#'                     "ENSMUSG00000030590", "ENSMUSG00000026188", "ENSMUSG00000014301",
#'                     "ENSMUSG00000073491", "ENSMUSG00000014529", "ENSMUSG00000036960",
#'                     "ENSMUSG00000058748", "ENSMUSG00000047388", "ENSMUSG00000002204",
#'                     "ENSMUSG00000034285", "ENSMUSG00000109129", "ENSMUSG00000035275",
#'                     "ENSMUSG00000051184", "ENSMUSG00000034424", "ENSMUSG00000041828",
#'                     "ENSMUSG00000029416", "ENSMUSG00000030468", "ENSMUSG00000029911",
#'                     "ENSMUSG00000055633", "ENSMUSG00000027495", "ENSMUSG00000029624",
#'                     "ENSMUSG00000045518", "ENSMUSG00000074259", "ENSMUSG00000035228",
#'                     "ENSMUSG00000038533", "ENSMUSG00000030401", "ENSMUSG00000014602",
#'                     "ENSMUSG00000041827", "ENSMUSG00000042345", "ENSMUSG00000028530",
#'                     "ENSMUSG00000038722", "ENSMUSG00000075088", "ENSMUSG00000039629",
#'                     "ENSMUSG00000067567", "ENSMUSG00000057594", "ENSMUSG00000005907",
#'                     "ENSMUSG00000027496")
#'   myTopAnatObject <- topAnat(myTopAnatData, geneList)
#'   resFis <- runTest(myTopAnatObject, algorithm = 'elim', statistic = 'fisher')
#'   ## Format results
#'   tableOver <- makeTable(myTopAnatData, myTopAnatObject, resFis, 0.1)
#' }
#' @import stats
#' @export

makeTable <- function(topAnatData, topAnatObject, results, cutoff=1){
  ## Perform some checks on the input data
  if(is.na(score(results)) || length(score(results)) == 0){
    stop("Problem: the results object is empty.")
  }
  if (!is.numeric(cutoff)){
    cutoff <- as.numeric(cutoff)
  }
  if(is.null(topAnatObject)){
    stop("Problem: the topAnatObject is empty.")
  }
  if( length(topAnatData$organ.names[,1]) == 0 ) {
    stop("Problem: the organ.names data frame of your topAnatData object is empty.")
  }

  ## retrieve p-values for the enrichment
  scores <- score(results)
  fdr <- p.adjust(p=scores, method = "fdr")
  topTerms <- sort(scores[fdr <= cutoff])
  topTerms <- as.data.frame(topTerms)

  if( nrow(topTerms) != 0 ){
    cat(paste0("\nBuilding the results table for the ", nrow(topTerms), " significant terms at FDR threshold of ", cutoff, "... "))
    odds <- termStat(topAnatObject, row.names(topTerms))
    foldEnrichment <- odds[2]/odds[3]

    ## Rounding odds, P-values and FDR to 2 decimal places
    ## Now commented because might be problematic for downstream reuse of the data
    # foldEnrichment <- format(foldEnrichment, digits=3)
    # topTerms <- format(topTerms, digits=3)
    # fdr[row.names(topTerms)] <- format(fdr[row.names(topTerms)], digits=3)

    topTerms <- cbind(odds, foldEnrichment, topTerms, fdr[row.names(topTerms)])
    names(topTerms) <- c("annotated", "significant", "expected", "foldEnrichment" , "pValue", "FDR")
    topTable <- merge(topAnatData$organ.names, topTerms, by.x=1, by.y=0)
    topTable <- topTable[order(as.numeric(topTable$pValue)), ]

    cat("Done\n\n")
    return(topTable)
  } else {
    cat("\nWarning: there was no significant term at a FDR threshold of", cutoff, "\n\n")
    return(NA)
  }
}
