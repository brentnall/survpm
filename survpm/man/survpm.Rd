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
This creates an object to calculate performance measures for prognostic survival models. In the base case the model may just be for the sub-hazard related to the endpoint of interest. The methods may also be applied to competing risks. For example, where there is competing mortality and a disease. Or more generally, with competing risks. It assumes discrete hazards from the risk model over time, but the survival time is continuous.
}
\usage{
survpm(df, m, maxT)
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
The maximum number of time periods.
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
mysumdta<-data.frame(myid=seq(1,21), myt=seq(0,20)+0.3, mycause=rep(c(0,1,2), 7), cbind(matrix(rep(seq(1,21)*0.03,each=21),ncol=21), matrix(rep(seq(1,11, by=0.5)*0.05,each=21),ncol=21)))
colnames(mysumdta)<-c("id", "t", "d", paste("H1-", 1:21, sep=""), paste("H2-", 1:21, sep=""))
## create spm object
myspm<-survpm(mysumdta, 2, 21)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}% use one of  RShowDoc("KEYWORDS")
