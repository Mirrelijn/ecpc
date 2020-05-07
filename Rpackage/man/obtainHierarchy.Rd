\name{obtainHierarchy}
\alias{obtainHierarchy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Obtain hierarchy
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function obtains the grouping on group level that defines the hierarchy;
if a group of covariates g is a subset of group h, then group h is an ancestor of group g (higher up in the hierarchy).
This hierarchy is used in adaptively discretising continuous co-data.
}
\usage{
obtainHierarchy(grouping, penalty = "LOG")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grouping}{
%%     ~~Describe \code{grouping} here~~
Grouping of groups of covariates with nested groups.
}
  \item{penalty}{
%%     ~~Describe \code{penalty} here~~
Default: "LOG" for a latent overlapping group approach (currently the only option in ecpc)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
We use the latent overlapping group (LOG) lasso penalty to define the hierarchical constraints as described in (Yan, Bien et al. 2007);
for each group g of covariates, we make a group on group level with group number g and the group numbers of its ancestors in the hierarchical tree. 
This way, group g can be selected if and only if all its ancestors are selected.
This function assumes that if group g is a subset of group h, then group h is an ancestor of group g. 
Note that this assumption does not necessarily hold for all hierarchies. The grouping on group level should then be coded manually.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A grouping on group level defining the hierarchy.
}
\references{
%% ~put references to the literature/web site here ~
Yan, X., Bien, J. et al. (2017). Hierarchical sparse modeling: A choice of two group lasso formulations. Statistical Science 32 531-560.

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
\code{\link{splitMedian}} to obtain a grouping of nested groups for continuous co-data.
}
\examples{
cont.codata <- seq(0,1,length.out=20) #continuous co-data
grouping <- splitMedian(values=cont.codata,split="lower",minGroupSize=5) #only split at lower continous co-data group
grouping.grouplvl <- obtainHierarchy(grouping) #obtain groups on group level defining the hierarchy

}