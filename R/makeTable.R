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
#' @param ordering A numeric indicating which column should be used to sort the data frame. If the column number is preceded by a \"-\" sign, results are displayed in decreasing ordering. Default is "7", returning data frame sorted by p-values in increasing order.
#'
#' @return A data frame with significantly enriched anatomical structures, sorted by p-value.
#'
#' @author Julien Roux
#'
#' @examples{
#' bgee <- Bgee$new(species="Bos_taurus", dataType="rna_seq")
#' myTopAnatData <- loadTopAnatData(bgee, stage="UBERON:0000092")
#' geneList <- as.factor(c(rep(0, times=85), rep(1, times=15)))
#' names(geneList) <- c("ENSBTAG00000000011","ENSBTAG00000000014","ENSBTAG00000000016",
#'                      "ENSBTAG00000000026","ENSBTAG00000000039","ENSBTAG00000000040",
#'                      "ENSBTAG00000000042","ENSBTAG00000000050","ENSBTAG00000000056",
#'                      "ENSBTAG00000000064","ENSBTAG00000000067","ENSBTAG00000000071",
#'                      "ENSBTAG00000000072","ENSBTAG00000000080","ENSBTAG00000000081",
#'                      "ENSBTAG00000000084","ENSBTAG00000000091","ENSBTAG00000000099",
#'                      "ENSBTAG00000000111","ENSBTAG00000000123","ENSBTAG00000000132",
#'                      "ENSBTAG00000000153","ENSBTAG00000000162","ENSBTAG00000000163",
#'                      "ENSBTAG00000000169","ENSBTAG00000000179","ENSBTAG00000000197",
#'                      "ENSBTAG00000000199","ENSBTAG00000000202","ENSBTAG00000000203",
#'                      "ENSBTAG00000000204","ENSBTAG00000000213","ENSBTAG00000000215",
#'                      "ENSBTAG00000000223","ENSBTAG00000000224","ENSBTAG00000000225",
#'                      "ENSBTAG00000000236","ENSBTAG00000000250","ENSBTAG00000000251",
#'                      "ENSBTAG00000000252","ENSBTAG00000000253","ENSBTAG00000000261",
#'                      "ENSBTAG00000000274","ENSBTAG00000000277","ENSBTAG00000000279",
#'                      "ENSBTAG00000000285","ENSBTAG00000000286","ENSBTAG00000000287",
#'                      "ENSBTAG00000000289","ENSBTAG00000000297","ENSBTAG00000000305",
#'                      "ENSBTAG00000000312","ENSBTAG00000000328","ENSBTAG00000000335",
#'                      "ENSBTAG00000000341","ENSBTAG00000000343","ENSBTAG00000000354",
#'                      "ENSBTAG00000000355","ENSBTAG00000000356","ENSBTAG00000000365",
#'                      "ENSBTAG00000000372","ENSBTAG00000000379","ENSBTAG00000000380",
#'                      "ENSBTAG00000000382","ENSBTAG00000000396","ENSBTAG00000000404",
#'                      "ENSBTAG00000000405","ENSBTAG00000000406","ENSBTAG00000000411",
#'                      "ENSBTAG00000000425","ENSBTAG00000000434","ENSBTAG00000000435",
#'                      "ENSBTAG00000000438","ENSBTAG00000000448","ENSBTAG00000000451",
#'                      "ENSBTAG00000000454","ENSBTAG00000000456","ENSBTAG00000000457",
#'                      "ENSBTAG00000000459","ENSBTAG00000000462","ENSBTAG00000000469",
#'                      "ENSBTAG00000000470","ENSBTAG00000000484","ENSBTAG00000000497",
#'                      "ENSBTAG00000000501","ENSBTAG00000009707","ENSBTAG00000026266",
#'                      "ENSBTAG00000021992","ENSBTAG00000005353","ENSBTAG00000005333",
#'                      "ENSBTAG00000006424","ENSBTAG00000026972","ENSBTAG00000010799",
#'                      "ENSBTAG00000010799","ENSBTAG00000014614","ENSBTAG00000014614",
#'                      "ENSBTAG00000045757","ENSBTAG00000046332","ENSBTAG00000046332",
#'                      "ENSBTAG00000008394")
#' myTopAnatObject <-  topAnat(myTopAnatData, geneList)
#' resFis <- runTest(myTopAnatObject, algorithm = 'elim', statistic = 'fisher')
#' ## Format results
#' tableOver <- makeTable(myTopAnatData, myTopAnatObject, resFis, 0.1)
#' }
#' @import stats
#' @export

makeTable <- function(topAnatData, topAnatObject, results, cutoff=1, ordering=7){
  ## Perform some checks on the input data
  if(any(is.na(score(results))) || length(score(results)) == 0){
    stop("ERROR: the results object is empty.")
  }
  if (!is.numeric(cutoff)){
    cutoff <- as.numeric(cutoff)
  }
  if(is.null(topAnatObject)){
    stop("ERROR: the topAnatObject is empty.")
  }
  if( length(topAnatData$organ.names[,1]) == 0 ) {
    stop("ERROR: the organ.names data frame of your topAnatData object is empty.")
  }

  dec <- FALSE
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }
  if(length(ordering) == 0){
    ordering <- 7
  }
  if(length(ordering) == 1 & is.numeric(ordering) & ordering <= 8){
    if (ordering < 0){
      dec <- TRUE
      ordering <- -ordering
    }
    if (!is.wholenumber(ordering)){
      stop("ERROR: problem with ordering argument. It should be an integer (column number).")
    }
  }
  else {
    stop("ERROR: problem with ordering argument. It should be an integer (column number), preceded by a \"-\" sign if decreasing ordering is needed. It should not be bigger than the number of columns in results data frame.")
  }
  ## retrieve p-values for the enrichment
  scores <- score(results)
  fdr <- p.adjust(p=scores, method = "fdr")
  topTerms <- sort(scores[fdr <= cutoff])
  topTerms <- as.data.frame(topTerms)

  if( nrow(topTerms) != 0 ){
    cat(paste0("\nBuilding the results table for the ", nrow(topTerms), " significant terms at FDR threshold of ", cutoff, "...\n"))
    odds <- termStat(topAnatObject, row.names(topTerms))
    foldEnrichment <- odds[2]/odds[3]

    ## Rounding odds, P-values and FDR to 2 decimal places
    ## Now commented because might be problematic for downstream reuse of the data. Maybe reintroduce with optional argument?
    # foldEnrichment <- format(foldEnrichment, digits=3)
    # topTerms <- format(topTerms, digits=3)
    # fdr[row.names(topTerms)] <- format(fdr[row.names(topTerms)], digits=3)

    topTerms <- cbind(odds, foldEnrichment, topTerms, fdr[row.names(topTerms)])
    topTable <- merge(topAnatData$organ.names, topTerms, by.x=1, by.y=0)
    names(topTable) <- c("organId", "organName", "annotated", "significant", "expected", "foldEnrichment" , "pValue", "FDR")

    cat(paste0("Ordering results by ", names(topTable)[ordering], " column in ", ifelse(dec, "decreasing", "increasing")," order...\n"))
    topTable <- topTable[order(topTable[, ordering], decreasing=dec), ]

    ## remove arbitrary numbers that are used as row names
    row.names(topTable) <- topTable$organId

    cat("Done\n\n")
    return(topTable)
  } else {
    cat("\nWARNING: there was no significant term at FDR threshold of", cutoff, "\n\n")
    return(NA)
  }
}
