\name{survpm}
\alias{survpm}
\alias{spm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~
Define a prognostic survival model object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This creates an object to calculate performance measures for prognostic survival models. The analysis may be used for survival analysis or competing risks sub-hazard models. Risk models are assumed to provide constant hazards over certain periods or time, or be reasonably approximated for performance analaysis when this is assumed. The survival times may be in continuous time. 
}
\usage{
survpm(df, m, maxT, ctime = FALSE, adjform = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
A data frame with the cohort survival and risk assessments.
}
  \item{m}{
%%     ~~Describe \code{m} here~~
The number of risks (e.g. survival analysis m=1).
}
  \item{maxT}{
%%     ~~Describe \code{maxT} here~~
The maximum number of time periods (model taken to have constant hazard in each time interval).
}
  \item{ctime}{
Whether to do analysis through follow-up time or not.
  }
  \item{adjform}{
Additional terms to add to Poisson regression for follow-up through time model (if wanted). 
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A survPM S4 class object.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Adam Brentnall
}
\note{
%%  ~~further notes~~
In development.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## demo data
mysumdta<-data.frame(myid=seq(1,21), myt=seq(0,20)+0.3, mycause=rep(c(0,1,2), 7), cbind(matrix(rep(seq(1,21)*0.03,each=21),ncol=21), matrix(rep(seq(1,11, by=0.5)*0.05,each=21),ncol=21)), grp1=c(rep(0,10), rep(1,11)), grp2=as.factor(rep(c(1,2,3), each=7)))

colnames(mysumdta)<-c("id", "t", "d", paste("H1-", 1:21, sep=""), paste("H2-", 1:21, sep=""), "x1", "x2")

myspm<-survpm(mysumdta, 2, 21)

## x1 + x2 adds this formula to the bespoke Poisson regression for time-dependent follow-up
myspm<-survpm(mysumdta, 2, 21, TRUE, "x1 + x2")
summary(myspm)

plot(myspm, 2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}% use one of  RShowDoc("KEYWORDS")
