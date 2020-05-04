\name{splitMedian}
\alias{splitMedian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Adaptively discretise continuous data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
splitMedian(values, index, depth, minGroupSize = 50, first = T, split = "left")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{values}{
%%     ~~Describe \code{values} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
}
  \item{depth}{
%%     ~~Describe \code{depth} here~~
}
  \item{minGroupSize}{
%%     ~~Describe \code{minGroupSize} here~~
}
  \item{first}{
%%     ~~Describe \code{first} here~~
}
  \item{split}{
%%     ~~Describe \code{split} here~~
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
function (values, index, depth, minGroupSize = 50, first = T, 
    split = "left") 
{
    if (missing(depth)) {
        depth <- floor(log2(length(values)/minGroupSize))
    }
    medianX <- median(values)
    indSmall <- values <= medianX
    indLarge <- values > medianX
    if (first) {
        groups <- list(index, index[indSmall], index[indLarge])
    }
    else {
        groups <- list(index[indSmall], index[indLarge])
    }
    if (depth > 1) {
        if (split == "both") {
            return(c(groups, splitMedian(values[indSmall], index[indSmall], 
                depth - 1, minGroupSize = minGroupSize, first = F, 
                split = split), splitMedian(values[indLarge], 
                index[indLarge], depth - 1, minGroupSize = minGroupSize, 
                first = F, split = split)))
        }
        else if (split == "left") {
            return(c(groups, splitMedian(values[indSmall], index[indSmall], 
                depth - 1, first = F, split = split)))
        }
        else {
            stop("split should be one of both or left")
        }
    }
    else {
        LargeEnough <- sapply(groups, function(x) {
            length(x) >= minGroupSize
        })
        if (all(LargeEnough)) {
            return(groups)
        }
        else {
            return(NULL)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
