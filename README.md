**BgeeDB: an R package for datasets retrieval from Bgee database**
------------------------------------------------------------------

BgeeDB is a collection of functions for everyday usage of Bgee database ( <http://bgee.org/> ) directly into R environment. In this pilot version of package, it is possible to retrieve the annotation of many experiments and libraries and download processed RNA-seq or Affymetrix data into R, according to Bgee pipeline. Currently, Bgee database supports gene expression data of 17 species.

> This package performs: listing annotation files and downloading gene expression data available in current version of Bgee database.

#### **Installation**

##### Install via install\_github()

``` r
# install dependecies
install.packages( c("data.table", "RCurl","devtools"))

# install the package
library(devtools)
install_github("wirawara/BgeeDB",build_vignettes=FALSE)
```

### **How to use BgeeDB package**

#### *Running example for Mus musculus*

##### List available species in current version of Bgee

``` r
library(BgeeDB)
listBgeeSpecies()
```

##### Choose the species and create a new Bgee species object

From the list above, please choose one species (eg. "Homo\_sapiens", "Mus\_musculus",...) and platform ("rna\_seq" or "affymetrix").

``` r
bgee <- Bgee$new(species = "Mus_musculus", datatype = "rna_seq")
```

##### Retrieve annotation for Mus musculus RNA-seq data

The **bgee$get\_annotation()** will output the list of experiments and libraries currently available in Bgee for RNA-seq of Mus musculus. The **bgee$get\_annotation()** loads the annotation in R, but also creates the Mus musculus folder in your current path, where it saves the downloaded annotation locally, so you can use the annotation for later as well. For loading the annotation separately, please use **fread("nameofyourannotation.tsv")**.

``` r
# the path where your folder with annotation will be saved. The folder is named after your chosen species.
getwd()
annotation_bgee_mouse <- bgee$get_annotation()
# access the experiments
annotation.experiments <- annotation_bgee_mouse[[1]]
# check the experiments annotation
head(annotation.experiments)

# access the libraries
annotation.libraries <- annotation_bgee_mouse[[2]]
# check the libraries annotation
head(annotation.libraries)
```

##### Download the processed Mus musculus RNA-seq data

The **bgee$get\_data()** will download RPKMs and counts for Mus musculus from all available experiments in Bgee database as a list form (see below). In case of downloaded data from all experiments for Mus musculus, **bgee$get\_data()** will save the downloaded data in your current folder as an **Bgee\_all\_experiment\_counts\_RPKM.rds** object, for later usage. You can load .rds again later using **readRDS("Bgee\_all\_experiment\_counts\_RPKM.rds")** from your current folder.

``` r
# the path where your data will be saved. 
getwd()

# download all RPKMs and counts for Mus musculus
data_bgee_mouse <- bgee$get_data()

# the number of experiments downloaded from Bgee
length(data_bgee_mouse)
# check the data
sapply(data_bgee_mouse, head)
# see your first experiment
data_bgee_experiment1 <- data_bgee_mouse[[1]]
```

Moreover, you can choose only one experiment to be downloaded from Bgee database, as in the example below. The data is then saved in .tsv file in your current folder.

``` r
# download RPKMs and counts only for GSE30617 for Mus musculus
data_bgee_mouse_gse30617 <- bgee$get_data(experiment.id = "GSE30617")
```
