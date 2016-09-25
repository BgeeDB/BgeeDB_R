---
output: html_document
---
# BgeeDB: an R package for datasets retrieval from Bgee database
##### Andrea Komljenović, Julien Roux, Marc Robinson-Rechavi, Frédéric Bastian

```BgeeDB``` is a collection of functions to import data from the Bgee database (<http://bgee.org/>) directly into R, and to facilitate downstream analyses, such as gene set enrichment test based on expression of genes in anatomical structures. Bgee provides annotated and processed expression data and expression calls from curated wild-type healthy samples, from humans and many animals.
 
The package retrieves the annotation of RNA-seq or Affymetrix experiments integrated into the Bgee database, and downloads into R the data reprocessed by the Bgee pipeline. It works for all the species in Bgee. The package also allows to run GO-like enrichment analyses based on anatomical terms, where genes are mapped to anatomical terms by expression patterns. This is the same as the TopAnat web-service available at (<http://bgee.org/?page=top_anat#/>), but with more flexibility in the choice of parameters and developmental stages, and is based on the ```topGO``` package.

This package allows: 
* 1. Listing annotation files gene of expression data available in the current version of Bgee database
* 2. Downloading the processed gene expression data available in the current version of Bgee database
* 3. Downloading the gene expression calls and annotations and using them to perform TopAnat analyses 

## Installation

### Install via Bioconductor

In R:
``` {r}
source("https://bioconductor.org/biocLite.R")
biocLite("BgeeDB")
```

*Warning:* you will be installing the package as it is in the `bioconductor` branch of the project, not the `master` branch. There could be minor differences, in particular you will need the very last version of R installed.
 
### Install via classic install

In the terminal:

    git clone https://github.com/wirawara/BgeeDB.git

Then in R:
``` {r}
install.packages("./BgeeDB", repos = NULL, type="source")
```
Or, download the project (`master` branch) by clicking the `Download ZIP` button on the web interface, and unzip it.

Then in R:
``` {r}
install.packages("./BgeeDB-master", repos = NULL, type="source")
```

### Install via install\_github()

``` {r}
# install the package
install.packages("devtools") # if you don't have devtools installed
library(devtools) 
install_github("wirawara/BgeeDB", build_vignettes=FALSE)
```


## How to use BgeeDB package

### Load the package
``` {r}
library(BgeeDB)
```

### Running example for Mus musculus

#### List available species in current version of Bgee
The ```listBgeeSpecies()``` function has several parameters. For example, it is possible to list species from different release versions of Bgee database with parameter ```release```, and order species according to offered columns with parameter ```ordering```.

``` {r}
listBgeeSpecies(release = "13.2", order = 2)
```

#### Choose the species and create a new Bgee species object

From the list above, please choose one species (e.g., "Homo\_sapiens", "Mus\_musculus",...) and platform ("rna\_seq" or "affymetrix"). The user can define the path where the data should be downloaded with parameter ```pathToData```, and the release version ```release```.

``` {r}
bgee <- Bgee$new(species = "Mus_musculus", datatype = "rna_seq", release = "13.2")
```

#### Retrieve annotation for Mus musculus RNA-seq data

The ```bgee$get_annotation()``` will output the list of experiments and libraries currently available in Bgee for RNA-seq of Mus musculus for the current version of Bgee database. The ```bgee$get_annotation()``` loads the annotation in R, but also creates the versioned Mus musculus folder in your current path, where it saves the downloaded annotation locally, so you can use the annotation for later as well.

``` {r}
# the path where your folder with annotation will be saved. The folder is named after your chosen species + version.
getwd()
annotation_bgee_mouse <- bgee$get_annotation()
# head the experiments and libraries
lapply(annotation_bgee_mouse, head)
```

#### Download the processed RNA-seq data for Mus musculus

The ```bgee$get_data()``` will download read counts and RPKMs for Mus musculus from all available experiments in Bgee database as a list (see below). In case of downloaded data from all experiments for Mus musculus, ```bgee$get_data()``` will save the downloaded data in your current folder for later usage. 

``` {r}
# the path where your data will be saved. 
getwd()

# download all RNAseq data for Mus musculus
data_bgee_mouse <- bgee$get_data()

# the number of experiments downloaded from Bgee
length(data_bgee_mouse)
# check the data
lapply(data_bgee_mouse, head)
# see your first experiment
data_bgee_experiment1 <- data_bgee_mouse[[1]]
```

*Note*: TPMs are going to be available in Bgee >=14.0. 


Alternatively, you can choose to download only one experiment from Bgee, as in the example below. The data is then saved as a .tsv file in your current folder.

``` {r}
# download RPKMs and counts only for GSE30617 for Mus musculus
data_bgee_mouse_gse30617 <- bgee$get_data(experiment.id = "GSE30617")
```

The function ```bgee$format_data()``` creates an ExpressionSet object including the expression data matrix, the annotation to Ensembl genes and the samples anatomical structure and stage annotation into (assayData, featureData and phenoData slots). This function also allows to filter out genes that are not called present in a given sample (giving them NA values). 

```{r}
# only present calls and rpkm values
gene.expression.mouse.rpkm <- bgee$format_data(data_bgee_mouse_gse30617, calltype = "present", stats = "rpkm")
gene.expression.mouse.rpkm 
```

#### Download the data to perform GO-like enrichment test for anatomical terms

The ```loadTopAnatData()``` function loads a mapping from genes to anatomical structures based on calls of expression in anatomical structures. It also loads the structure of the anatomical ontology and the names of anatomical structures.

```{r}
# Loading calls of expression based on RNA-seq data only
myTopAnatData <- loadTopAnatData(species="Mus_musculus", datatype="rna_seq")
# Loading calls observed in embryonic stages only
myTopAnatData <- loadTopAnatData(species="Mus_musculus", datatype="rna_seq", stage="UBERON:0000068")

# Look at the data
lapply(myTopAnatData, head)
```

*Note*: the results are stored in files (see the ```pathToData``` arguments). To save time, if you query again with the exact same parameters, these cache files will be read instead of querying the web-service. So do not delete the files in the working folder if you plan to perform additional queries.

#### Prepare a topGO object allowing to perform GO-like enrichment test for anatomical terms, for Mus musculus

First we need to prepare a list of genes in the foreground and in the background. The input format is the same as the gene list required to build a ```topGOdata``` object in the ```topGO``` package: a vector with background genes as names, and 0 or 1 values depending if a gene is in the foreground or not. In this example we will look at genes, annotated with "spermatogenesis" in the Gene Ontology (using the ```biomaRt``` package). 

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)

# Foreground genes are those with a GO annotation "spermatogenesis"
myGenes <- getBM(attributes= "ensembl_gene_id", filters=c("go_id"), values=list(c("GO:0007283")), mart=ensembl)

# Background are all genes with a GO annotation
universe <- getBM(attributes= "ensembl_gene_id", filters=c("with_go_go"), values=list(c(TRUE)), mart=ensembl)

# Prepares the gene list vector 
geneList <- factor(as.integer(universe[,1] %in% myGenes[,1]))
names(geneList) <- universe[,1]
head(geneList)
summary(geneList == 1)

# Prepares the topGO object
myTopAnatObject <-  topAnat(myTopAnatData, geneList)
```

*Warning*: This can be long, especially if the gene list is large, since the anatomical ontology is large and expression calls will be propagated through the whole ontology (e.g., expression in the forebrain will also be counted as expression in parent structures such as the brain, nervous system, etc). Consider running a script in batch mode if you have multiple analyses to do.

#### Launch an enrichment test for anatomical terms

For this step, see the vignette of the ```topGO``` package for more details, as you have to directly use the tests implemented in the ```topGO``` package, as shown in this example:

```{r}
results <- runTest(myTopAnatObject, algorithm = 'classic', statistic = 'fisher')
# You can also use the topGO decorrelation methods, for example the "weight" method to get less redundant results
results <- runTest(myTopAnatObject, algorithm = 'weight', statistic = 'fisher')
```

*Warning*: This can be long because of the size of the ontology. Consider running a script in batch mode if you have multiple analyses to do.

#### Format the table of results after an enrichment test for anatomical terms

We built the ```makeTable``` function to filter and format the test results. Results are sorted by p-value, and FDR values are calculated. 

```{r}
# Display results sigificant at a 1% FDR threshold
tableOver <- makeTable(myTopAnatData, myTopAnatObject, results, cutoff = 0.01)
# Display all results
tableOver <- makeTable(myTopAnatData, myTopAnatObject, results, cutoff = 1)
```

*Warning*: it is debated whether FDR correction is appropriate on enrichment test results, since tests on different terms of the ontologies are not independent. A nice discussion can be found in the vignette of the ```topGO``` package.

