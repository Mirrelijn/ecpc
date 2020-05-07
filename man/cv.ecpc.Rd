\name{cv.ecpc}
\alias{cv.ecpc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Cross-validation for ecpc
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Cross-validates ecpc and returns model fit, summary statistics and cross-validated performance measures.
}
\usage{
cv.ecpc(Y,X,postselection="elnet+dense",maxsel=30,
                    type.measure="MSE",outerfolds=10,
                    ncores=1,balance=T,lambdas,model,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
Response data; n-dimensional vector (n: number of samples) for linear and logistic outcomes, or \code{\link{Surv}} object for Cox survival.
}
  \item{X}{
%%     ~~Describe \code{X} here~~
Observed data; (nxp)-dimensional matrix (p: number of covariates) with each row the observed high-dimensional feature vector of a sample.
}
  \item{postselection}{
%%     ~~Describe \code{postselection} here~~
Type of posterior selection method used to obtain a parsimonious model of maxsel covariates, or FALSE if no parsimonious model is needed. Possible options are "elnet,dense" (default), "elnet,sparse", "BRmarginal,dense", "BRmarginal,sparse" or "DSS".
}
  \item{maxsel}{
%%     ~~Describe \code{maxsel} here~~
Maximum number of covariates to be selected a posteriori, in addition to all unpenalised covariates. If maxsel is a vector, multiple parsimonious models are returned.
}
  \item{type.measure}{
  Type of cross-validated performance measure returned.
}
  \item{outerfolds}{
  Number of cross-validation folds.
}
  \item{ncores}{
  Number of cores used for parallel processing the outer cross-validation folds.
}
  \item{balance}{
  (logistic, Cox) Should folds be balanced in response?
}
  \item{lambdas}{
  A vector of global ridge penalties for each fold may be given.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
Type of model for the response; linear, logistic or Cox.
}
  \item{...}{
  Additional arguments used in \code{\link{ecpc}}.
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
A list with the following elements:
\item{ecpc.fit}{List with the ecpc model fit in each fold}
\item{dfPred}{Data frame with information about out-of-bag predictions}
\item{dfGrps}{Data frame with information about estimated group and grouping weights across folds}
\item{dfCVM}{Data frame with cross-validated performance metric}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


}
