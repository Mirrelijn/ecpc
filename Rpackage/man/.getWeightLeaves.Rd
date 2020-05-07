\name{.getWeightLeaves}
\alias{.getWeightLeaves}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Find local variance of leaf groups
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns a data frame with the local variance in the leaf groups in a hierarchical tree,
instead of the group weights per group.
}
\usage{
.getWeightLeaves(dfGrps, Grouping, grouping.grouplvl, values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfGrps}{
%%     ~~Describe \code{dfGrps} here~~
Dataframe with the following variables:
  -Group: factor with group names
  -Group.weight: group weight of each group
  -Fold: number indicating which fold in the cross-validation is used
}
  \item{Grouping}{
%%     ~~Describe \code{Grouping} here~~
List containing covariate indices for each group in the grouping.
}
  \item{grouping.grouplvl}{
%%     ~~Describe \code{grouping.grouplvl} here~~
(optional) Groups of groups, e.g. in a hierarchical structure. 
}
  \item{values}{
%%     ~~Describe \code{values} here~~
(optional) Values of continuous co-data.
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
df: dataframe similar to 'dfGrps' but with group weights for leaf groups
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