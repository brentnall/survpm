\name{SurvPM-class}
\Rdversion{1.1}
\docType{class}
\alias{SurvPM-class}
\alias{calcH,SurvPM-method}
\alias{calcHt,SurvPM-method}
\alias{calcOE,SurvPM-method}
\alias{initialize,SurvPM-method}

\title{Class \code{"SurvPM"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
The class related to performance analysis of a prognostic risk model. The model is specified such that the cumulative sub-hazards of different competing risks over discrete time periods may be obtained (e.g. yearly), while survival to the first of (potentially) several competing risks is continuous. It keeps the source data and derived data needed for the analysis, in addition to methods to undertake the analysis.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SurvPM", crData, m, maxT)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{crData}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{crhaz}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{crH}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{m}:}{Object of class \code{"integer"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{maxT}:}{Object of class \code{"integer"} ~~ }
    \item{\code{pmH}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{resOEtot}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{crHt}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{calcH}{\code{signature(object = "SurvPM")}: ... }
    \item{calcHt}{\code{signature(object = "SurvPM")}: ... }
    \item{calcOE}{\code{signature(object = "SurvPM")}: ... }
    \item{initialize}{\code{signature(.Object = "SurvPM")}: ... }
    \item{plot}{\code{signature(x = "SurvPM")}: ... }
    \item{summary}{\code{signature(object = "SurvPM")}: ... }
    \item{tdepcal}{\code{signature(object = "SurvPM")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
Adam Brentnall
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("SurvPM")
}
\keyword{classes}
