\name{NEWS}
\title{News for package \pkg{BgeeDB}}

\section{Version 2.5.0 (2018-02)}{
  \itemize{
    \item Compatibility with new tar.gz compressed download files.
    \item Better management of potential issues with the ontology structure.
    \item Better management of .RData object datasets used to create the vignette.
    \item Fixed small issue with species names composed of more than 2 names (e.g Canis lupus familiaris)
  }
}

\section{Version 2.3.2 (2017-05)}{
  \itemize{
    \item Compatibility of the TopAnat part of the package with the new Bgee release 14 (www.bgee.org/bgee14/).
    \item Better management of potential issues with the ontology structure
  }
}

\section{Version 2.3.1 (2017-05)}{
  \itemize{
    \item Compatibility with the new Bgee release 14 (www.bgee.org/bgee14/) only for the processed data download part.
    \item The TopAnat part still temporarily uses data from Bgee release 13.2, until the webservice is ready.
    \item Retro compatibility with Bgee release 13.2 (www.bgee.org/bgee13/)
  }
}

\section{Version 2.2.0 (2017-04)}{
  \itemize{
    \item Instead of returning a topGO-compatible object, topAnat.R now returns an object from the topAnatData class, an extension of topGOdata class from the topGO package.
    \item Fixed small issue with management of data types given as input by the user (dataType argument when creating new Bgee class)
    \item Fixed bug in experiment Id check step. Now accomodates SRA Ids.
    \item Fixed data frames header names that included double dots.
    \item Removed dependency to biomaRt in vignette. Code is still detailed but not run, instead a pre-created gene list object is loaded from the data/ directory.
  }
}

\section{Version 2.0.0 (2016-10-10)}{
  \itemize{
    \item Implemented possibility to deal with different Bgee releases.
    \item Improved storage and versioning of cached files.
    \item Implemented use of API key to query our servers in order to prevent overloading and spamming.
    \item Improved management of downloading errors.
    \item Harmonized the use of a Bgee class object by all functions of the package. For example, loadTopAnatData() now requires an input Bgee class object to specify species, dataType and pathToData arguments.
    \item Added input Bgee class object to output of loadtopAnatData() function.
    \item Created new getAnnotation(), getData() and formatData() independent functions to replace the Bgee class methods get_annotation(), get_data() and format_data().
    \item In formatData() function, when affymetrix data is used, the "stats" parameter is automatically set to "intensities".
    \item Added possibility to reproduce an analysis offline if all data files were previously downloaded in cache.
    \item Fixed data frames headers including spaces to more convenient headers with spaces replaced by dots.
    \item Harmonized use of camelCase in functions arguments.
    \item Added argument allowing to sort result table in makeTable() function.
    \item Implemented management of TPMs as expression unit in future Bgee releases.
    \item Updated vignette.
  }
}

\section{Version 1.0.3 (2016-08-31)}{
  \itemize{
    \item Update of format_data() function to output an Expression Set object.
    \item Fixed makeTable "FDR" column which was a factor instead of a numeric.
    \item Fixed get_data() and format_data() functions, which did not work when multiple chip types were available for an experiment.
  }
}

\section{Version 1.0.0 (2016-05-03)}{
  \itemize{
    \item Public release of BgeeDB Bioconductor package. Enjoy!
  }
}
