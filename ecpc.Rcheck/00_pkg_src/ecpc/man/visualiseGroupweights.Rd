\name{visualiseGroupweights}
\alias{visualiseGroupweights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Visualise estimated group weights
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plot group weights from multiple cross-validation folds.
}
\usage{
visualiseGroupweights(dfGrps, Grouping, grouping.grouplvl, values, 
                      widthBoxplot = 0.05, boxplot = T, jitter = T, 
                      ps = 1.5, ls = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfGrps}{
%%     ~~Describe \code{dfGrps} here~~
Data frame containing the following variables;
  'Group': factor with group names;
  'Group.weight': group weight of each group;
  'Fold': number indicating which fold in the cross-validation is used.
}
  \item{Grouping}{
%%     ~~Describe \code{Grouping} here~~
List of G elements containing covariate indices for each group
}
  \item{grouping.grouplvl}{
%%     ~~Describe \code{grouping.grouplvl} here~~
(optional): groups on group level, e.g. defining a hierarchical structure. 
}
  \item{values}{
%%     ~~Describe \code{values} here~~
(optional): values of continuous co-data. If given, group weights are plotted against these value.
}
  \item{widthBoxplot}{
%%     ~~Describe \code{widthBoxplot} here~~
Width of boxplot.
}
  \item{boxplot}{
%%     ~~Describe \code{boxplot} here~~
Should a boxplot be plotted?
}
  \item{jitter}{
%%     ~~Describe \code{jitter} here~~
Should point estimates be jittered?
}
  \item{ps}{
%%     ~~Describe \code{ps} here~~
Point size.
}
  \item{ls}{
%%     ~~Describe \code{ls} here~~
Line size.
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
\code{\link{visualiseGrouping}} to visualise groupings and \code{\link{visualiseGroupingweights}} to plot estimated grouping weights.
}
\examples{
# #required libraries
# library(ggplot2)
# library(dplyr)
#
# #discrete groups
# grouping1 <- list(1:20,21:40)
# dfGrps1 <- data.frame(Group=as.factor(rep(c(1,2),each=10)),
#                      Group.weight=c(rnorm(10,0.5,0.01),rnorm(10,2,0.05)),
#                      Fold=rep(1:10,2))
# visualiseGroupweights(dfGrps1, Grouping=grouping1)
# 
# #continous co-data groups
# cont.codata <- seq(0,1,length.out=40) #continuous co-data
# grouping2 <- splitMedian(values=cont.codata,split="lower",minGroupSize=10) #only split at lower continous co-data group
# grouping.grouplvl <- obtainHierarchy(grouping2) #obtain groups on group level defining the hierarchy
# 
# #simulate random group weights around 1
# dfGrps2 <- data.frame(Group=as.factor(rep(1:length(grouping2),each=10)),
#                       Group.weight=c(rnorm(10*length(grouping2),1,0.01)),
#                       Fold=rep(1:10,length(grouping2))) 
# visualiseGroupweights(dfGrps2, Grouping=grouping2, grouping.grouplvl=grouping.grouplvl) #plot group weights per group
# visualiseGroupweights(dfGrps2, Grouping=grouping2, grouping.grouplvl=grouping.grouplvl, 
#                       values=cont.codata) #plot group weights per leaf group in the hierarchical tree
}
