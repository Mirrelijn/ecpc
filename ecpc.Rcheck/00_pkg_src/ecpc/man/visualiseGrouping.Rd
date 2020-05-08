\name{visualiseGrouping}
\alias{visualiseGrouping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Visualise a grouping
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Visualises a grouping in a graph, with directed edges indicating the hierarchy.
}
\usage{
visualiseGrouping(Grouping, groupweights, grouping.grouplvl, nodeSize = 10, ls = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Grouping}{
%%     ~~Describe \code{Grouping} here~~
List of G groups of covariates.
}
  \item{groupweights}{
%%     ~~Describe \code{groupweights} here~~
(optional) vector with G group weights; if given, group weights are visualised too.
}
  \item{grouping.grouplvl}{
%%     ~~Describe \code{grouping.grouplvl} here~~
List of G_2 groups defining a hierarchy.
}
  \item{nodeSize}{
%%     ~~Describe \code{nodeSize} here~~
Size of the nodes in the visualisation; scalar.
}
  \item{ls}{
%%     ~~Describe \code{ls} here~~
Line size; scalar.
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
A ggplot object.

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
\code{\link{visualiseGroupingweights}} to plot estimated grouping weights. and \code{\link{visualiseGroupweights}} to plot estimated group weights.
}
\examples{
# #required libraries
# library(Matrix)
# library(dplyr)
# library(igraph)
# library(ggplot2)
# library(ggraph)

# #groups without hierarchical constraints
# grouping <- list("Group1"=c(1:20),"Group2"=c(15,30))
# visualiseGrouping(grouping,c(0.5,2))

# #hierarchical groups
# cont.codata <- seq(0,1,length.out=20) #continuous co-data
# hierarchicalgrouping <- splitMedian(values=cont.codata,split="lower",minGroupSize=5) #only split at lower continous co-data group
# grouping.grouplvl <- obtainHierarchy(hierarchicalgrouping) #obtain groups on group level defining the hierarchy
# visualiseGrouping(hierarchicalgrouping, grouping.grouplvl=grouping.grouplvl)

}
