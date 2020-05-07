\name{visualiseGroupingweights}
\alias{visualiseGroupingweights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Visualise estimated grouping weights
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plot grouping weights from multiple cross-validation folds.
}
\usage{
visualiseGroupingweights(dfGrps, GroupingNames, hist = F, boxplot = T, 
                          jitter = T, ps = 1.5, width = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfGrps}{
%%     ~~Describe \code{dfGrps} here~~
Data frame containing the following variables;
  'Grouping': factor with grouping names;
  'Grouping.weight': groupings weight of each grouping;
  'Fold': number indicating which fold in the cross-validation is used.
}
  \item{GroupingNames}{
%%     ~~Describe \code{GroupingNames} here~~
Vector with names of the groupings.
}
  \item{hist}{
%%     ~~Describe \code{hist} here~~
Should histogram be plotted?
}
  \item{boxplot}{
%%     ~~Describe \code{boxplot} here~~
Should boxplot be used or points?
}
  \item{jitter}{
%%     ~~Describe \code{jitter} here~~
Should grouping weights be jittered?
}
  \item{ps}{
%%     ~~Describe \code{ps} here~~
Point size.
}
  \item{width}{
%%     ~~Describe \code{width} here~~
Width of jitter.
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
Plot in ggplot object.
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
\code{\link{visualiseGrouping}} to visualise groupings and \code{\link{visualiseGroupweights}} to plot estimated group weights.
}
\examples{
# #required libraries
# library(ggplot2)
# library(dplyr)
#
# dfGrps <- data.frame(Grouping=rep(c(1,2),each=10),
#                      Grouping.weight=c(rnorm(10,0,0.01),rnorm(10,1,0.05)),
#                      Fold=rep(1:10,2))
# GroupingNames <- c("Grouping1","Grouping2")
# visualiseGroupingweights(dfGrps, GroupingNames, hist = F, boxplot = T,jitter=T)
}

