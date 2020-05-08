\name{postSelect}
\alias{postSelect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Perform posterior selection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given data and estimated parameters from a previously fit multi-group ridge penalised model, perform posterior selection to find a parsimonious model.
}
\usage{
postSelect(X, Y, beta, intrcpt = 0, penfctr, postselection = c("elnet+dense",
           "elnet+sparse","BRmarginal,dense","BRmarginal,sparse","DSS"), 
           maxsel = 30, penalties, model=c("linear","logistic","cox"), 
           tauglobal, sigmahat = 1, muhatp = 0, X2 = NaN, Y2 = NaN)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
Observed data: data of p penalised and unpenalised covariates on n samples; (nxp)-dimensional matrix.
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
Response data; n-dimensional vector (linear, logistic) or \code{\link{Surv}} object (Cox survival).
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
Estimated regression coefficients from the previously fit model.
}
  \item{intrcpt}{
%%     ~~Describe \code{intrcpt} here~~
Estimated intercept from the previously fit model.
}
  \item{penfctr}{
%%     ~~Describe \code{penfctr} here~~
As in glmnet penalty.factor; p-dimensional vector with a 0 if covariate is not penalised, 1 if covariate is penalised.
}
  \item{postselection}{
%%     ~~Describe \code{postselection} here~~
Posterior selection method to be used.
}
  \item{maxsel}{
%%     ~~Describe \code{maxsel} here~~
Maximum number of covariates to be selected a posteriori, in addition to all unpenalised covariates. If maxsel is a vector, multiple parsimonious models are returned.
}
  \item{penalties}{
%%     ~~Describe \code{penalties} here~~
Estimated multi-group ridge penalties for all penalised covariates from the previously fit model; vector of length the number of penalised covariates.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
Type of model for the response.
}
  \item{tauglobal}{
%%     ~~Describe \code{tauglobal} here~~
Estimated global prior variance from the previously fit model.
}
  \item{sigmahat}{
%%     ~~Describe \code{sigmahat} here~~
(linear model only) estimated variance parameter from the previously fit model.
}
  \item{muhatp}{
%%     ~~Describe \code{muhatp} here~~
(optional) Estimated multi-group prior means for the penalised covariates from the previously fit model.
}
  \item{X2}{
%%     ~~Describe \code{X2} here~~
(optional) Independent observed data.
}
  \item{Y2}{
%%     ~~Describe \code{Y2} here~~
(optional) Independent response data.
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
  \item{betaPost}{Estimated regression coefficients for parsimonious models. If 'maxsel' is a vector,
  'betaPost' is a matrix with each column the vector estimate corresponding to the maximum number of 
  selected covariates given in 'maxsel'.}
  \item{a0}{Estimated intercept coefficient for parsimonious models.}
  \item{YpredPost}{If independent test set 'X2' is given, posterior selection model predictions for the test set.}
  \item{MSEPost}{If independent test set 'X2', 'Y2' is given, mean squared error of the posterior selection model predictions.}
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