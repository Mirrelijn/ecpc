\name{simDat}
\alias{simDat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simDat(n, p, n2 = 20, muGrp, varGrp, indT, sigma = 1, model = "linear", flag = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{n2}{
%%     ~~Describe \code{n2} here~~
}
  \item{muGrp}{
%%     ~~Describe \code{muGrp} here~~
}
  \item{varGrp}{
%%     ~~Describe \code{varGrp} here~~
}
  \item{indT}{
%%     ~~Describe \code{indT} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{flag}{
%%     ~~Describe \code{flag} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, p, n2 = 20, muGrp, varGrp, indT, sigma = 1, model = "linear", 
    flag = F) 
{
    X <- mvrnorm(n, rep(0, p), diag(rep(1, p)))
    X2 <- mvrnorm(n2, rep(0, p), diag(rep(1, p)))
    meansX <- apply(X, 2, mean)
    Xctd <- t(t(X) - meansX)
    meansX2 <- apply(X2, 2, mean)
    X2ctd <- t(t(X2) - meansX2)
    beta <- mvrnorm(1, muGrp[indT], diag(varGrp[indT]))
    if (model == "linear") {
        Y <- Xctd \%*\% beta + rnorm(n, 0, sd = sigma)
        Y2 <- X2ctd \%*\% beta + rnorm(n2, 0, sd = sigma)
        if (flag) {
            print(paste("Simulated data for", model, "regression"))
            plot(Y)
            points(Xctd \%*\% beta, col = "red")
            legend("right", c("Y", "X*beta"), pch = c(1, 1), 
                col = c("black", "red"))
            hist(Y)
        }
    }
    if (model == "logistic") {
        expXb <- exp(Xctd \%*\% beta)
        Y <- rbinom(n, 1, expXb/(1 + expXb))
        expX2b <- exp(X2ctd \%*\% beta)
        Y2 <- rbinom(n2, 1, expX2b/(1 + expX2b))
        if (flag) {
            print(paste("Simulated data for", model, "regression"))
            plot(Y)
            points(expXb/(1 + expXb), col = "red")
            legend("right", c("Y", "expit(X*beta)"), pch = c(1, 
                1), col = c("black", "red"))
            hist(expXb/(1 + expXb))
        }
    }
    return(list(beta = beta, Xctd = Xctd, Y = Y, X2ctd = X2ctd, 
        Y2 = Y2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
