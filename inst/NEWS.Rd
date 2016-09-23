\name{NEWS}
\title{News for Package \pkg{BgeeDB}}

\section{Changes in version 1.0.4 (2016-09-24)}{
  \itemize{
    \item Implemented a way to deal with different Bgee releases
    \item Implemented a way to deal with presence of TPMs in future Bgee releases
    \item Harmonized the use of the species parameter between the Bgee object creation and the loadTopAnatData() function.
  }
}

\section{Changes in version 1.0.3 (2016-08-31)}{
  \itemize{
    \item Update of format_data to output an Expression Set object.
    \item Fixed makeTable "FDR" column which was a factor instead of a numeric.
    \item Fixed get_data() function bug when multiple chip types were available for an experiment.
  }
}


\section{Changes in version 0.99.6 (2016-04-28)}{
  \itemize{
    \item  Bug fix in format_data function.
  }
}

\section{Changes in version 0.99.4 (2016-04-26)}{
  \itemize{
    \item  Bug fix for Affymetrix data.
  }
}

\section{Changes in version 0.99.3 (2016-04-22)}{
  \itemize{
    \item  Added download for Affymetrix data.
  }
}


\section{Changes in version 0.0.4 (2016-03-20)}{
  \itemize{
    \item  Added loadTopAnatData function.
  }
}

\section{Changes in version 0.0.1 (2016-03-13)}{
  \itemize{
    \item Added topAnat script.
      This function produces a topAnatObject,
      ready to use for gene set enrichment testing
      using functions from the topGO package
  }
}
