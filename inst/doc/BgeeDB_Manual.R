## -----------------------------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly=TRUE))
    #install.packages("BiocManager")
#BiocManager::install("BgeeDB")

## ---- message = FALSE, warning = FALSE----------------------------------------
library(BgeeDB)

## -----------------------------------------------------------------------------
listBgeeSpecies()

## -----------------------------------------------------------------------------
listBgeeSpecies(release = "13.2", order = 2)

## -----------------------------------------------------------------------------
bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")

## -----------------------------------------------------------------------------
annotation_bgee_mouse <- getAnnotation(bgee)
# list the first experiments and libraries
lapply(annotation_bgee_mouse, head)

## -----------------------------------------------------------------------------
# download all RNA-seq experiments from mouse
data_bgee_mouse <- getData(bgee)
# number of experiments downloaded
length(data_bgee_mouse)
# check the downloaded data
lapply(data_bgee_mouse, head)
# isolate the first experiment
data_bgee_experiment1 <- data_bgee_mouse[[1]]

## -----------------------------------------------------------------------------
# download data for GSE30617
data_bgee_mouse_gse30617 <- getData(bgee, experimentId = "GSE30617")

## ----eval=FALSE---------------------------------------------------------------
#  # Examples of data downloading using different filtering combination
#  # retrieve mouse RNA-Seq data for heart or brain
#  data_bgee_mouse_filters <- getData(bgee, anatEntityId = c("UBERON:0000955","UBERON:0000948"))
#  # retrieve mouse RNA-Seq data for heart (UBERON:0000955) or brain (UBERON:0000948) part of the experiment GSE30617
#  data_bgee_mouse_filters <- getData(bgee, experimentId = "GSE30617", anatEntityId = c("UBERON:0000955","UBERON:0000948"))
#  # retrieve mouse RNA-Seq data for heart (UBERON:0000955) or brain (UBERON:0000948) from post-embryonic stage (UBERON:0000092)
#  data_bgee_mouse_filters <- getData(bgee, stageId = "UBERON:0000092", anatEntityId = c("UBERON:0000955","UBERON:0000948"))

## -----------------------------------------------------------------------------
# use only present calls and fill expression matric with FPKM values
gene.expression.mouse.fpkm <- formatData(bgee, data_bgee_mouse_gse30617, callType = "present", stats = "fpkm")
gene.expression.mouse.fpkm 

## -----------------------------------------------------------------------------
# Creating new Bgee class object
bgee <- Bgee$new(species = "Danio_rerio")

## -----------------------------------------------------------------------------
# Loading calls of expression
myTopAnatData <- loadTopAnatData(bgee)
# Look at the data
## str(myTopAnatData)

## ---- eval=FALSE--------------------------------------------------------------
#  ## Loading silver and gold expression calls from affymetrix data made on embryonic samples only
#  ## This is just given as an example, but is not run in this vignette because only few data are returned
#  bgee <- Bgee$new(species = "Danio_rerio", dataType="affymetrix")
#  myTopAnatData <- loadTopAnatData(bgee, stage="UBERON:0000068", confidence="silver")

## ---- eval=FALSE--------------------------------------------------------------
#  # if (!requireNamespace("BiocManager", quietly=TRUE))
#      # install.packages("BiocManager")
#  # BiocManager::install("biomaRt")
#  library(biomaRt)
#  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="drerio_gene_ensembl", host="mar2016.archive.ensembl.org")
#  
#  # get the mapping of Ensembl genes to phenotypes. It will corresponds to the background genes
#  universe <- getBM(filters=c("phenotype_source"), value=c("ZFIN"), attributes=c("ensembl_gene_id","phenotype_description"), mart=ensembl)
#  
#  # select phenotypes related to pectoral fin
#  phenotypes <- grep("pectoral fin", unique(universe$phenotype_description), value=T)
#  
#  # Foreground genes are those with an annotated phenotype related to "pectoral fin"
#  myGenes <- unique(universe$ensembl_gene_id[universe$phenotype_description %in% phenotypes])
#  
#  # Prepare the gene list vector
#  geneList <- factor(as.integer(unique(universe$ensembl_gene_id) %in% myGenes))
#  names(geneList) <- unique(universe$ensembl_gene_id)
#  summary(geneList)
#  
#  # Prepare the topGO object
#  myTopAnatObject <-  topAnat(myTopAnatData, geneList)

## -----------------------------------------------------------------------------
data(geneList)
myTopAnatObject <-  topAnat(myTopAnatData, geneList)

## -----------------------------------------------------------------------------
results <- runTest(myTopAnatObject, algorithm = 'weight', statistic = 'fisher')

## -----------------------------------------------------------------------------
# Display results sigificant at a 1% FDR threshold
tableOver <- makeTable(myTopAnatData, myTopAnatObject, results, cutoff = 0.01)
head(tableOver)

## -----------------------------------------------------------------------------
# In order to retrieve significant genes mapped to the term " paired limb/fin bud"
term <- "UBERON:0004357"
termStat(myTopAnatObject, term) 

# 198 genes mapped to this term for Bgee 14.0 and Ensembl 84
genesInTerm(myTopAnatObject, term)
# 48 significant genes mapped to this term for Bgee 14.0 and Ensembl 84
annotated <- genesInTerm(myTopAnatObject, term)[["UBERON:0004357"]]
annotated[annotated %in% sigGenes(myTopAnatObject)]

## ----eval = FALSE-------------------------------------------------------------
#  bgee <- Bgee$new(species="Mus_musculus", release = "14.1")
#  # delete all old .rds files of species Mus musculus
#  deleteOldData(bgee)

## ----eval = FALSE-------------------------------------------------------------
#  bgee <- Bgee$new(species="Mus_musculus", release = "14.1")
#  # delete local SQLite database of species Mus musculus from Bgee 14.1
#  deleteLocalData(bgee)
