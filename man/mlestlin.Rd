\name{mlestlin}
\alias{mlestlin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum marginal likelihood estimation for the linear model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Finds the maximum marginal likelihood estimates for the overall ridge penalty (prior variance) and noise level.
}
\usage{
mlestlin(Y, X, lambda = NaN, sigmasq = NaN, mu = NaN, tausq = NaN)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{sigmasq}{
%%     ~~Describe \code{sigmasq} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{tausq}{
%%     ~~Describe \code{tausq} here~~
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
function (Y, X, lambda = NaN, sigmasq = NaN, mu = NaN, tausq = NaN) 
{
    maxv <- var(Y)
    p <- dim(X)[2]
    n <- dim(X)[1]
    prms <- paste(c("l", "s", "m", "t")[is.nan(c(lambda, sigmasq, 
        mu, tausq))], collapse = "")
    if (prms == "t" || prms == "mt") {
        tausq <- sigmasq/lambda
        prms <- gsub("t", "", prms)
    }
    if (prms == "l" || prms == "lm") {
        lambda <- sigmasq/tausq
        prms <- gsub("l", "", prms)
    }
    if (prms == "s" || prms == "sm") {
        sigmasq <- lambda * tausq
        prms <- gsub("s", "", prms)
    }
    switch(prms, lsmt = {
        sim2 = function(ts) {
            tausq <- ts[1]
            sigmasq <- ts[2]
            mu <- ts[3]
            n <- nrow(X)
            varY <- X \%*\% t(X) * exp(tausq) + diag(rep(1, n)) * 
                exp(sigmasq)
            meanY <- X \%*\% rep(mu, p)
            mlk <- -dmvnorm(c(Y), mean = meanY, sigma = varY, 
                log = T)
            return(mlk)
        }
        op <- optim(c(log(0.01), log(maxv), 0), sim2)
        tausq <- exp(op$par[1])
        sigmasq <- exp(op$par[2])
        mu <- op$par[3]
        lambda <- sigmasq/tausq
    }, lsm = {
        sim2 = function(ts) {
            sigmasq <- ts[1]
            mu <- ts[2]
            n <- nrow(X)
            varY <- X \%*\% t(X) * exp(tausq) + diag(rep(1, n)) * 
                exp(sigmasq)
            meanY <- X \%*\% rep(mu, p)
            mlk <- -dmvnorm(c(Y), mean = meanY, sigma = varY, 
                log = T)
            return(mlk)
        }
        op <- optim(c(log(maxv), 0), sim2)
        sigmasq <- exp(op$par[1])
        mu <- op$par[2]
        lambda <- sigmasq/tausq
    }, lst = {
        sim2 = function(ts) {
            tausq <- ts[1]
            sigmasq <- ts[2]
            n <- nrow(X)
            varY <- X \%*\% t(X) * exp(tausq) + diag(rep(1, n)) * 
                exp(sigmasq)
            meanY <- X \%*\% rep(mu, p)
            mlk <- -dmvnorm(c(Y), mean = meanY, sigma = varY, 
                log = T)
            return(mlk)
        }
        op <- optim(c(log(0.01), log(maxv)), sim2)
        tausq <- exp(op$par[1])
        sigmasq <- exp(op$par[2])
        lambda <- sigmasq/tausq
    }, lmt = {
        sim2 = function(ts) {
            tausq <- ts[1]
            mu <- ts[2]
            n <- nrow(X)
            varY <- X \%*\% t(X) * exp(tausq) + diag(rep(1, n)) * 
                sigmasq
            meanY <- X \%*\% rep(mu, p)
            mlk <- -dmvnorm(c(Y), mean = meanY, sigma = varY, 
                log = T)
            return(mlk)
        }
        op <- optim(c(log(0.01), 0), sim2)
        tausq <- exp(op$par[1])
        mu <- op$par[2]
        lambda <- sigmasq/tausq
    }, lt = {
        sim2 = function(ts) {
            tausq <- ts[1]
            n <- nrow(X)
            varY <- X \%*\% t(X) * exp(tausq) + diag(rep(1, n)) * 
                sigmasq
            meanY <- X \%*\% rep(mu, p)
            mlk <- -dmvnorm(c(Y), mean = meanY, sigma = varY, 
                log = T)
            return(mlk)
        }
        op <- optim(c(log(0.01)), sim2, method = "Brent", lower = log(1e-10), 
            upper = log(100))
        tausq <- exp(op$par[1])
        lambda <- sigmasq/tausq
    }, smt = {
        sim2 = function(ts) {
            tausq <- log(exp(ts[1])/lambda)
            sigmasq <- ts[1]
            mu <- ts[2]
            n <- nrow(X)
            varY <- X \%*\% t(X) * exp(tausq) + diag(rep(1, n)) * 
                exp(sigmasq)
            meanY <- X \%*\% rep(mu, p)
            mlk <- -dmvnorm(c(Y), mean = meanY, sigma = varY, 
                log = T)
            return(mlk)
        }
        op <- optim(c(log(maxv), 0), sim2, method = "L-BFGS-B", 
            lower = c(log(1e-10), -10^{
                99
            }), upper = c(log(2 * maxv), 10^{
                99
            }))
        tausq <- exp(op$par[1])/lambda
        sigmasq <- exp(op$par[1])
        mu <- op$par[2]
    }, st = {
        sim2 = function(ts) {
            tausq <- log(exp(ts)/lambda)
            sigmasq <- ts
            n <- nrow(X)
            varY <- X \%*\% t(X) * exp(tausq) + diag(rep(1, n)) * 
                exp(sigmasq)
            meanY <- X \%*\% rep(mu, p)
            mlk <- -dmvnorm(c(Y), mean = meanY, sigma = varY, 
                log = T)
            return(mlk)
        }
        op <- optim(c(log(maxv)), sim2, method = "Brent", lower = log(1e-10), 
            upper = log(2 * maxv))
        tausq <- exp(op$par[1])/lambda
        sigmasq <- exp(op$par[1])
    }, ls = {
        sim2 = function(ts) {
            sigmasq <- ts[1]
            n <- nrow(X)
            varY <- X \%*\% t(X) * exp(tausq) + diag(rep(1, n)) * 
                sigmasq
            meanY <- X \%*\% rep(mu, p)
            mlk <- -dmvnorm(c(Y), mean = meanY, sigma = varY, 
                log = T)
            return(mlk)
        }
        op <- optim(c(log(maxv)), sim2, method = "Brent", lower = log(1e-10), 
            upper = log(2 * maxv))
        sigmasq <- exp(op$par[1])
        lambda <- sigmasq/tausq
    }, m = {
        sim2 = function(ts) {
            mu <- ts
            n <- nrow(X)
            varY <- X \%*\% t(X) * tausq + diag(rep(1, n)) * sigmasq
            meanY <- X \%*\% rep(mu, p)
            mlk <- -dmvnorm(c(Y), mean = meanY, sigma = varY, 
                log = T)
            return(mlk)
        }
        op <- optim(0, sim2, method = "Brent", lower = -10^{
            99
        }, upper = 10^{
            99
        })
        mu <- op$par[1]
    })
    mlests <- c(lambda, sigmasq, mu, tausq)
    return(mlests)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
