\name{survpm-package}
\alias{survpm-package}
\alias{calcH}
\alias{calcHt}
\alias{calcOE}
\alias{calcOEt}
\alias{fn.format}
\alias{fn.formatCI}
\alias{fn.formatp}
\alias{fn.spmsim}
\alias{checkSurvPM}
\docType{package}
\title{
\packageTitle{survpm}
}
\description{
\packageDescription{survpm}
}
\details{
\packageDESCRIPTION{survpm}
\packageIndices{survpm}
An overview of how to use the package, including the most important functions
}
\author{
\packageAuthor{survpm}

Maintainer: \packageMaintainer{survpm}
}
\references{
Literature or other references for background information 
}
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
maxT <- 20
haz1 <-  (1:20)/100
haz2 <-  rep(2,20)/100
myn<-200

mysumdta <- fn.spmsim(haz1, haz2, myn)
myspm<-survpm(mysumdta, m=2, maxT, ctime=TRUE)

par(mfrow=c(2,2))
plot(myspm); plot(myspm,2)
plot(myspm,1,2); plot(myspm,2,2)
summary(myspm)

}

