\name{produceFolds}
\alias{produceFolds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Produce folds
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Produce folds for cross-validation.
}
\usage{
produceFolds(nsam, outerfold, response, model = c("logistic","survival","other"), 
balance = TRUE, fixedfolds = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsam}{
%%     ~~Describe \code{nsam} here~~
Number of samples
}
  \item{outerfold}{
%%     ~~Describe \code{outerfold} here~~
Number of folds.
}
  \item{response}{
%%     ~~Describe \code{response} here~~
Response data.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
Type of model for the response.
}
  \item{balance}{
%%     ~~Describe \code{balance} here~~
Should folds be balanced in response?
}
  \item{fixedfolds}{
%%     ~~Describe \code{fixedfolds} here~~
Should fixed folds be used?
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
