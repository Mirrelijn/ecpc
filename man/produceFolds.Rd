\name{produceFolds}
\alias{produceFolds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Produce folds
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
produceFolds(nsam, outerfold, response, model = "logistic", balance = TRUE, fixedfolds = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsam}{
%%     ~~Describe \code{nsam} here~~
}
  \item{outerfold}{
%%     ~~Describe \code{outerfold} here~~
}
  \item{response}{
%%     ~~Describe \code{response} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{balance}{
%%     ~~Describe \code{balance} here~~
}
  \item{fixedfolds}{
%%     ~~Describe \code{fixedfolds} here~~
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
function (nsam, outerfold, response, model = "logistic", balance = TRUE, 
    fixedfolds = F) 
{
    if (fixedfolds) 
        set.seed(3648310)
    if (!balance) {
        rand <- sample(1:nsam)
        grs1 <- floor(nsam/outerfold)
        grs2 <- grs1 + 1
        ngr1 <- outerfold * grs2 - nsam
        folds <- lapply(1:outerfold, function(xg) {
            if (xg <= ngr1) 
                els <- rand[(1 + (xg - 1) * grs1):(xg * grs1)]
            else els <- rand[(ngr1 * grs1 + 1 + (xg - ngr1 - 
                1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
            return(els)
        })
    }
    else {
        if (model == "logistic") 
            if (class(response) == "factor") 
                nev <- which((as.numeric(response) - 1) == 1)
            else nev <- which(response == 1)
        if (model == "survival") 
            nev <- which(response[, 1] == 1)
        nsamev <- length(nev)
        randev <- sample(nev)
        grs1 <- floor(nsamev/outerfold)
        grs2 <- grs1 + 1
        ngr1 <- outerfold * grs2 - nsamev
        foldsev <- lapply(1:outerfold, function(xg) {
            if (xg <= ngr1) 
                els <- randev[(1 + (xg - 1) * grs1):(xg * grs1)]
            else els <- randev[(ngr1 * grs1 + 1 + (xg - ngr1 - 
                1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
            return(els)
        })
        nonev <- setdiff(1:nsam, nev)
        nsamnonev <- length(nonev)
        randnonev <- sample(nonev)
        grs1 <- floor(nsamnonev/outerfold)
        grs2 <- grs1 + 1
        ngr1 <- outerfold * grs2 - nsamnonev
        foldsnonev <- lapply(1:outerfold, function(xg) {
            if (xg <= ngr1) 
                els <- randnonev[(1 + (xg - 1) * grs1):(xg * 
                  grs1)]
            else els <- randnonev[(ngr1 * grs1 + 1 + (xg - ngr1 - 
                1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
            return(els)
        })
        folds <- lapply(1:outerfold, function(i) c(foldsev[[i]], 
            foldsnonev[[i]]))
    }
    return(folds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
