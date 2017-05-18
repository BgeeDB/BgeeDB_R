#' @title Produces an object allowing to perform GO-like enrichment of anatomical terms using the topGO package
#'
#' @description This function produces a topAnatObject, ready to use for gene set enrichment testing using functions from the topGO package. This object uses the Uberon ontology instead of the GO ontology.
#'
#' @details To perform the enrichment test for expression in anatomical structures for each term of Uberon ontology (browsable at \url{http://www.ontobee.org/ontology/UBERON}), the data are formatted to use the topGO package for testing. This package is interesting because it propagates the mapping of gene to terms to parent terms, and it possesses a pannel of enrichment tests and decorrelation methods. Expert users should be able to use information from the topAnatObject to test enrichment with other packages than topGO.
#'
#' @param topAnatData a list including a gene2anatomy list, an organ.relationships list and an organ.names data.frame, produced by the function loadTopAnatData().
#'
#' @param geneList Vector indicating foreground and background genes. Names of the vector indicate the background genes. Values are 1 (gene in foreground) or 0 (gene not in foreground).
#'
#' @param nodeSize Minimum number of genes mapped to a node for it to be tested. Default is 10.
#'
#' @param ... Additional parameters as passed to build topGOdata object in topGO package.
#'
#' @return topAnatObject, a topAnatData class object, ready for gene set enrichment testing with topGO.
#'
#' @author Julien Roux
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'   myTopAnatData <- loadTopAnatData(bgee)
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
#'   myTopAnatObject <- topAnat(myTopAnatData, geneList, nodeSize=1)
#' }
#'
#' @import topGO graph
#' @export

topAnat <- function(topAnatData, geneList, nodeSize = 10, ... ){

  ## Test if topAnatData not empty
  cat("\nChecking topAnatData object.............\n")
  if( length(topAnatData$gene2anatomy) == 0 ) {
    stop("ERROR: the gene2anatomy list of your topAnatData object is empty.")
  }
  if( length(topAnatData$organ.relationships) == 0 ) {
    stop("ERROR: the organ.relationships list of your topAnatData object is empty.")
  }
  if( length(topAnatData$organ.names[,1]) == 0 ) {
    stop("ERROR: the organ.names data frame of your topAnatData object is empty.")
  }

  ## Test if gene list is fine
  cat("\nChecking gene list......................\n")
  if (!is.factor(geneList)){
    geneList <- as.factor(geneList)
  }
  if (length(geneList) == 0){
    stop("ERROR: the gene list provided is empty.")
  }
  if (length(levels(geneList)) != 2){
    stop("ERROR: the gene list provided is not in the right format (should be a named vector including only 0 and 1 values).")
  }
  if (length(geneList) < 100) {
    cat("\nWARNING: Given the low number of genes provided, it is very unlikely that the test will have enough power.\n")
  }
  ## If geneList includes genes not present in topAnatData$gene2anatomy, restrict to these genes
  if (sum(names(geneList) %in% names(topAnatData$gene2anatomy)) != length(geneList)){
    cat("\nWARNING: Some genes in your gene list have no expression data in Bgee, and will not be included in the analysis.", sum(names(geneList) %in% names(topAnatData$gene2anatomy)), "genes in background will be kept.\n")
  }

  if (nodeSize == 0){
    stop("ERROR: the node size parameter has to be at least 1.")
  }
  if (!is.numeric(nodeSize)){
    stop("ERROR: the node size parameter should be a numeric, integer value.")
  }
  nodeSize <- as.integer(nodeSize)

  ## Building the modified topAnatData object. This also reports to the user how many genes are in the background / foreground
  topAnatObject <- new("topAnatData",
                         description = "topAnatData class object, ready for anatomical ontology enrichment test",
                         ontology = "UBERON ontology describing animal anatomical structures",
                         allGenes = geneList,
                         nodeSize = nodeSize,
                         parentMapping = topAnatData$organ.relationships,
                         gene2Nodes = topAnatData$gene2anatomy
                         )

  ## Create a hash from the topAnatObject and the topAnatData objects
  myAnalysisInfo <- c(topAnatObject, topAnatData)
  ## Use library digest to create a SHA512 hash, that will be used as analysis Id
  analysisId <- digest(myAnalysisInfo, algo = "sha512")
  ## cat(paste0("\nAnalysis Id hash built: ", analysisId, "\n"))

  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ## The analysis id is used for our internal statistics. It is a secure hash that does not allow the
  ## identification of the gene lists and datasets used by the user, but will inform us on the
  ## number of distinct analyses that are run with the package.
  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ## Send a query to our webservice for statistics purposes of topAnat use (even when cached files are used)
  if (topAnatData$bgee.object$sendStats == TRUE){
    myUrl <- paste0(topAnatData$bgee.object$topAnatUrl, "?page=stats&action=launch_top_anat_analysis&api_key=", topAnatData$bgee.object$apiKey, "&analysis_id=", analysisId, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    ## Send query, but no need to wait for the answer
    try(getURL(myUrl, followLocation = TRUE, .opts = list(timeout = 1)), silent=TRUE)
  }

  return(topAnatObject)
}

#################################################################
## Sets of functions written with help of Adrian Alexa (author of the topGO package, personnal communication) to make topGO work with another ontology than the Gene Ontology

.annFUN.gene2Nodes <- function(feasibleGenes = NULL, gene2Nodes) {
  ## Restrict the mappings to the feasibleGenes set
  if(!is.null(feasibleGenes))
    gene2Nodes <- gene2Nodes[intersect(names(gene2Nodes), feasibleGenes)]

  ## Throw-up the genes which have no annotation
  if(any(is.na(gene2Nodes)))
    gene2Nodes <- gene2Nodes[!is.na(gene2Nodes)]

  gene2Nodes <- gene2Nodes[sapply(gene2Nodes, length) > 0]

  ## Get all the Terms and keep a one-to-one mapping with the genes
  allNodes <- unlist(gene2Nodes, use.names = FALSE)
  geneID <- rep(names(gene2Nodes), sapply(gene2Nodes, length))

  return(split(geneID, allNodes))
}
#################################################################

.buildGraph.topology <- function(knownNodes, parentMapping) {
  ## first build the lookUp table for the terms
  nodeLookUp <- new.env(hash = TRUE, parent = emptyenv())

  ## warping functions for a easier acces to the lookUp table
  isNodeInDAG <- function(node) {
    return(exists(node, envir = nodeLookUp, mode = 'logical', inherits = FALSE))
  }
  setNodeInDAG <- function(node) {
    assign(node, TRUE, envir = nodeLookUp)
  }

  ## get the root of the ontology
  GENE.ONTO.ROOT <- setdiff(unique(unlist(parentMapping, use.names = FALSE)),
      names(parentMapping))
  if (length(GENE.ONTO.ROOT) > 1){
    cat('\n Multiple roots were found... Possibly some nodes are not listed in the organ.names data frame\n')
  }

  ## we read all the database once and the access the list
  adjLookUP <- as.list(parentMapping)

  ## we use an environment of environments to store edges: (this way is faster)
  ## in the end we will coerce it to a list of list and build a graphNEL obj.
  edgeEnv <- new.env(hash = TRUE, parent = emptyenv())

  ## add the arc (u --> v) to edgeEnv of type :
  envAddEdge <- function(u, v) {
    assign(v, 0, envir = get(u, envir = edgeEnv))
  }

  ## recursivly build the induced graph starting from one node
  buildInducedGraph <- function(node) {
    ## if we have visited the node, there is nothing to do
    if(isNodeInDAG(node))
      return(1)

    ## we put the node in the graph and we get his parents
    setNodeInDAG(node)    # we visit the node
    assign(node, new.env(hash = TRUE, parent = emptyenv()), envir = edgeEnv) # adj list

    if(node == GENE.ONTO.ROOT)
      return(2)

    adjNodes <- adjLookUP[[node]]

    ## debuging point! should not happen!
    if(length(adjNodes) == 0)
      cat('\n There are no adj nodes for node: ', node, '\n')

    for(i in 1:length(adjNodes)) {
      x <- as.character(adjNodes[i])
      envAddEdge(node, x)
      buildInducedGraph(x)
    }

    return(0)
  }

  ## we start from the most specific nodes
  lapply(knownNodes, buildInducedGraph)

  ## now we must transform our env into a Graph structure
  ## for now we use lapply, later we can do it with eapply
  .graphNodes <- ls(edgeEnv)
  .edgeList <- eapply(edgeEnv,
                      function(adjEnv) {
                        aux <- as.list(adjEnv)
                        return(list(edges = match(names(aux), .graphNodes),
                                    weights = as.numeric(aux)))
                      })

  ## now we can build the graphNEL object
  return(new('graphNEL',
             nodes = .graphNodes,
             edgeL = .edgeList,
             edgemode = 'directed'))
}
######################## topAnatData class ########################
## This class is an extension of the topGOdata class from
## the topGO package. The additional fields are:
## - parentMapping: a list describing parent-child
##   relationships in the ontology to be used in place of
##   the Gene Ontology
## - gene2Nodes: a list providing the mapping of genes to
##   term from the new ontology

## Declare topAnatData as an extension of topGOdata class
setClass("topAnatData",
         representation = representation(
           ## Added slot to give us the structure of the ontology
           parentMapping = "character",
           ## Added slot to give us the mapping of genes to terms of the ontology
           gene2Nodes = "character"),
         contains = "topGOdata"
)

setMethod("initialize", "topAnatData",
          function(.Object,
                   ## which Ontology to be used
                   ontology = "Uberon",
                   ## a named numeric or factor, the names are the genes ID
                   allGenes,
                   ## function to select the signif. genes
                   geneSelectionFun = NULL,
                   ## description of the class
                   description = character(0),
                   ## minimum node size
                   nodeSize = 1,
                   ## the child-parent relationship
                   parentMapping,
                   ## annotation data
                   gene2Nodes,
                   ## additional parameters
                   ...) {

            .Object@description <- description
            .Object@ontology <- ontology

            ## some checking
            if(is.null(names(allGenes)))
              stop("allGenes must be a named vector")

            if(!is.factor(allGenes) && !is.numeric(allGenes))
              stop("allGenes should be a factor or a numeric vector")

            .Object@allGenes <- names(allGenes)

            if(is.factor(allGenes)) {
              if(length(levels(allGenes)) != 2)
                stop("allGenes must be a factor with 2 levels")
              .Object@allScores <- factor(as.character(allGenes))
              .Object@geneSelectionFun <- function(x) {
                return(as.logical(as.integer(levels(x)))[x])
              }
            }
            else {
              .Object@allScores <- as.numeric(allGenes)

              ## function to select which genes are significant
              if(is.null(geneSelectionFun))
                warning("No function to select the significant genes provided!")
              .Object@geneSelectionFun <- geneSelectionFun
            }

            ## size of the nodes which will be pruned
            .Object@nodeSize = as.integer(max(nodeSize, 1))

            ## this function is returning a list of terms from the specified ontology
            ## with each entry being a vector of genes
            cat("\nBuilding most specific Ontology terms... ")
            mostSpecificTerms <- .annFUN.gene2Nodes(feasibleGenes = .Object@allGenes, gene2Nodes = gene2Nodes)
            cat(" ( ", length(mostSpecificTerms), " Ontology terms found. )\n")

            ## the graph is build starting from the most specific terms
            cat("\nBuilding DAG topology................... ")
            g <- .buildGraph.topology(names(mostSpecificTerms), parentMapping)
            cat(" ( ",  graph::numNodes(g), " Ontology terms and ", graph::numEdges(g), " relations. )\n")

            ## probably is good to store the leves but for the moment we don't
            .nodeLevel <- buildLevels(g, leafs2root = TRUE)

            ## annotate the nodes in the graph with genes
            cat("\nAnnotating nodes (Can be long).......... ")
            g <- topGO:::mapGenes2GOgraph(g, mostSpecificTerms, nodeLevel = .nodeLevel) ## leafs2root

            ## select the feasible genes
            gRoot <- getGraphRoot(g)
            feasibleGenes <- ls(graph::nodeData(g, n = gRoot, attr = "genes")[[gRoot]])
            cat(" ( ", length(feasibleGenes), " genes annotated to the Ontology terms. )\n")

            .Object@feasible <- .Object@allGenes %in% feasibleGenes


            ## prune the GO graph
            if(.Object@nodeSize > 1) {
              cc <- .countsInNode(g, graph::nodes(g))
              .Object@graph <- graph::subGraph(names(cc)[cc >= .Object@nodeSize], g)
            } else {
              .Object@graph <-  g
            }

            .Object
          })
