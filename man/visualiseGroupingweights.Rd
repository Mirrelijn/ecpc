\name{visualiseGroupingweights}
\alias{visualiseGroupingweights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Visualise estimated grouping weights
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
visualiseGroupingweights(dfGrps, GroupingNames, hist = F, boxplot = T, jitter = T, ps = 1.5, width = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfGrps}{
%%     ~~Describe \code{dfGrps} here~~
}
  \item{GroupingNames}{
%%     ~~Describe \code{GroupingNames} here~~
}
  \item{hist}{
%%     ~~Describe \code{hist} here~~
}
  \item{boxplot}{
%%     ~~Describe \code{boxplot} here~~
}
  \item{jitter}{
%%     ~~Describe \code{jitter} here~~
}
  \item{ps}{
%%     ~~Describe \code{ps} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
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
function (dfGrps, GroupingNames, hist = F, boxplot = T, jitter = T, 
    ps = 1.5, width = 0.5) 
{
    if (!missing(GroupingNames) && length(GroupingNames) > 0) {
        if (is.factor(dfGrps$Grouping)) {
            dfGrps$Grouping <- factor(dfGrps$Grouping, levels = levels(dfGrps$Grouping)[levels(dfGrps$Grouping) \%in\% 
                unique(dfGrps$Grouping)], labels = GroupingNames)
        }
        else {
            dfGrps$Grouping <- factor(GroupingNames[dfGrps$Grouping], 
                levels = GroupingNames, labels = GroupingNames)
        }
    }
    else {
        dfGrps$Grouping <- as.factor(dfGrps$Grouping)
    }
    if (hist) {
        p1 <- ggplot(dfGrps) + aes(x = Grouping.weight) + geom_histogram(aes(fill = Grouping), 
            position = "identity", bins = 20, alpha = 0.6) + 
            scale_fill_discrete(name = "Grouping") + labs(x = "Prior variance grouping weight") + 
            theme_bw() + theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14), 
            axis.title.y = element_text(size = 14), legend.text = element_text(size = 12), 
            legend.title = element_text(size = 14))
    }
    else if (!(boxplot | jitter)) {
        SummdfGrps <- dfGrps \%>\% group_by(Grouping) \%>\% summarise(meanWeight = mean(Grouping.weight), 
            q25Weight = quantile(Grouping.weight, 0.25), q75Weight = quantile(Grouping.weight, 
                0.75), maxWeight = max(Grouping.weight), minWeight = min(Grouping.weight)) \%>\% 
            ungroup()
        p1 <- ggplot(SummdfGrps) + geom_linerange(aes(x = Grouping, 
            ymax = q75Weight, ymin = q25Weight, col = Grouping), 
            linetype = "dashed") + geom_point(aes(x = Grouping, 
            y = minWeight, col = Grouping), shape = 2) + geom_point(aes(x = Grouping, 
            y = meanWeight, col = Grouping), shape = 1) + geom_point(aes(x = Grouping, 
            y = maxWeight, col = Grouping), shape = 6) + scale_color_discrete(name = "Grouping") + 
            labs(y = "Prior variance grouping weight", x = "Grouping") + 
            theme_bw() + theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14), 
            axis.title.y = element_text(size = 14), legend.text = element_text(size = 12), 
            legend.title = element_text(size = 14))
    }
    else {
        dfGrps2 <- dfGrps \%>\% group_by(Fold) \%>\% distinct(Grouping, 
            .keep_all = T) \%>\% ungroup()
        p1 <- ggplot(dfGrps2) + scale_color_discrete(name = "Grouping") + 
            labs(y = "Prior variance grouping weight", x = "Grouping") + 
            theme_bw() + theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14), 
            axis.title.y = element_text(size = 14), legend.text = element_text(size = 12), 
            legend.title = element_text(size = 14))
        if (boxplot) {
            p1 <- p1 + geom_boxplot(aes(x = Grouping, y = Grouping.weight, 
                col = Grouping), outlier.shape = NA, width = width)
        }
        if (jitter) {
            p1 <- p1 + geom_jitter(aes(x = Grouping, col = Grouping, 
                y = Grouping.weight), size = ps, alpha = 0.6, 
                height = 0, width = width/2)
        }
    }
    return(p1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
