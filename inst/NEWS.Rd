\name{NEWS}
\title{News for package \pkg{BgeeDB}}

\section{Version 1.2.0 (2016-10-10)}{
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
