\name{hierarchicalLasso}
\alias{hierarchicalLasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fit hierarchical lasso using LOG penalty
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Fits a linear regression model penalised with a hierarchical lasso penalty, using a latent overlapping group (LOG) lasso penalty.
}
\usage{
hierarchicalLasso(X, Y, grouping)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
nxp matrix with observed data
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
nx1 vector with response data
}
  \item{grouping}{
%%     ~~Describe \code{grouping} here~~
list with hierarchical group indices
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The LOG penalty can be used to impose hierarchical constraints in the estimation of regression coefficients (Yan, Bien et al. 2007), e.g. a group of covariates (child node in the hierarchical tree) may be selected only if another group is selected (parent node in the hierarchical tree).
This function uses the simple implementation for the LOG penalty described in (Jacob, Obozinski and Vert, 2009). Faster and more scalable algorithms may be available but not yet used in this pacakage.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list with the following elements;
  \item{betas}{Estimated regression coefficients.}
  \item{a0}{Estimated intercept.}
  \item{lambda}{Estimated penalty parameter.}
  \item{group.weights}{Fixed group weights used in the LOG-penalty.}
}
\references{
%% ~put references to the literature/web site here ~
Yan, X., Bien, J. et al. (2017). Hierarchical sparse modeling: A choice of two group lasso formulations. Statistical Science 32 531–560.

Jacob, L., Obozinski, G. and Vert, J.-P. (2009). Group lasso with overlap and graph lasso. In: Proceedings of the 26th annual international conference on machine learning 433-440. ACM.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
