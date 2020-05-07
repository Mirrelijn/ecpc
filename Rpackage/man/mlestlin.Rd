\name{mlestlin}
\alias{mlestlin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum marginal likelihood estimation for the linear model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Finds the maximum marginal likelihood estimates for the global ridge penalty (prior variance) and noise level for the linear model.
}
\usage{
mlestlin(Y, X, lambda = NaN, sigmasq = NaN, mu = NaN, tausq = NaN)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
Linear response.
}
  \item{X}{
%%     ~~Describe \code{X} here~~
Observed data.
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
Ridge penalty; if given, lamba is fixed.
}
  \item{sigmasq}{
%%     ~~Describe \code{sigmasq} here~~
Variance parameter; if given, sigmasq is fixed.
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
Prior mean/target for regression coefficients; if given fixed (default set to 0).
}
  \item{tausq}{
%%     ~~Describe \code{tausq} here~~
Prior variance parameter; equals sigmasq/lambda.
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
A vector with the maximum marginal likelihood estimates for 'lambda', 'sigmasq', 'mu' and 'tausq'.
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