\name{NEWS}
\title{News for Package \pkg{BgeeDB}}

\section{Changes in version 1.2 (2016-10-10)}{
  \itemize{
    \item Implemented a way to deal with different Bgee releases and versioning of cached files.
    \item Implemented a way to deal with addition of TPMs as expression unit future Bgee releases.
    \item Fixed data frames headers including spaces to more convenient headers with spaces replaced by dots.
    \item Harmonized the use of a Bgee class object by other functions throughout the packages.
    \itam Created new getAnnotation(), getData() and formatData() independent functions to replace the Bgee class methods get_annotation(), get_data() and format_data().
    \item In formatData() function, when affymetrix data is used, the "stats" parameter is automatically set to "intensities".
    \item Added argument allowing to sort result table in makeTable() function
  }
}

\section{Changes in version 1.0.3 (2016-08-31)}{
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
