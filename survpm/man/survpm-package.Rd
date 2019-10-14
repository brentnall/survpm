\name{survpm-package}
\alias{survpm-package}
\alias{calcH}
\alias{calcHt}
\alias{calcOE}
\alias{calcOEt}
\alias{fn.format}
\alias{fn.formatCI}
\alias{fn.formatp}
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
## demo data
mysumdta<-data.frame(myid=seq(1,21), myt=seq(0,20)+0.3, mycause=rep(c(0,1,2), 7), cbind(matrix(rep(seq(1,21)*0.03,each=21),ncol=21), matrix(rep(seq(1,11, by=0.5)*0.05,each=21),ncol=21)))
colnames(mysumdta)<-c("id", "t", "d", paste("H1-", 1:21, sep=""), paste("H2-", 1:21, sep=""))
## create spm object
myspm<-survpm(mysumdta, 2, 21)
## plot through time
plot(myspm, idx=2)
}
