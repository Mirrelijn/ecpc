\name{visualiseGroupweights}
\alias{visualiseGroupweights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Visualise estimated group weights
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
visualiseGroupweights(dfGrps, Grouping, grouping.grouplvl, values, widthBoxplot = 0.05, boxplot = T, jitter = T, ps = 1.5, ls = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfGrps}{
%%     ~~Describe \code{dfGrps} here~~
}
  \item{Grouping}{
%%     ~~Describe \code{Grouping} here~~
}
  \item{grouping.grouplvl}{
%%     ~~Describe \code{grouping.grouplvl} here~~
}
  \item{values}{
%%     ~~Describe \code{values} here~~
}
  \item{widthBoxplot}{
%%     ~~Describe \code{widthBoxplot} here~~
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
  \item{ls}{
%%     ~~Describe \code{ls} here~~
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
function (dfGrps, Grouping, grouping.grouplvl, values, widthBoxplot = 0.05, 
    boxplot = T, jitter = T, ps = 1.5, ls = 1) 
{
    if (!missing(Grouping) && length(names(Grouping)) > 0) {
        if (is.factor(dfGrps$Group)) {
            dfGrps$Group <- factor(dfGrps$Group, levels = levels(dfGrps$Group)[levels(dfGrps$Group) \%in\% 
                unique(dfGrps$Group)], labels = names(Grouping))
        }
        else {
            dfGrps$Group <- factor(names(Grouping)[dfGrps$Group])
        }
    }
    if (missing(values) || length(values) == 0) {
        SummdfGrps <- dfGrps \%>\% group_by(Group) \%>\% summarise(meanWeight = mean(Group.weight), 
            q25Weight = quantile(Group.weight, 0.25), q75Weight = quantile(Group.weight, 
                0.75), maxWeight = max(Group.weight), minWeight = min(Group.weight)) \%>\% 
            ungroup()
        if (!(boxplot | jitter)) {
            p1 <- ggplot(SummdfGrps) + geom_hline(yintercept = 1, 
                col = "grey", linetype = "dashed", size = ls, 
                alpha = 0.6) + geom_linerange(aes(x = Group, 
                ymax = q75Weight, ymin = q25Weight, col = Group), 
                linetype = "dashed") + geom_point(aes(x = Group, 
                y = minWeight, col = Group), shape = 2) + geom_point(aes(x = Group, 
                y = meanWeight, col = Group), shape = 1) + geom_point(aes(x = Group, 
                y = maxWeight, col = Group), shape = 6) + labs(y = "Prior variance weight") + 
                theme_bw() + theme(axis.text.x = element_text(size = 12), 
                axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14), 
                axis.title.y = element_text(size = 14), legend.text = element_text(size = 12), 
                legend.title = element_text(size = 14))
        }
        else {
            p1 <- ggplot(dfGrps) + geom_hline(yintercept = 1, 
                col = "grey", linetype = "dashed", size = ls, 
                alpha = 0.6) + labs(y = "Prior variance weight") + 
                theme_bw() + theme(axis.text.x = element_text(size = 12), 
                axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14), 
                axis.title.y = element_text(size = 14), legend.text = element_text(size = 12), 
                legend.title = element_text(size = 14))
            if (boxplot) {
                p1 <- p1 + geom_boxplot(aes(x = Group, y = Group.weight, 
                  col = Group), outlier.shape = NA, show.legend = F)
            }
            if (jitter) {
                p1 <- p1 + geom_jitter(aes(x = Group, y = Group.weight, 
                  col = Group), show.legend = F, size = ps, alpha = 0.6, 
                  height = 0, width = widthBoxplot/2)
            }
        }
        return(p1)
    }
    else {
        getIndBetas <- function(p, ind) {
            indBetas <- lapply(1:p, function(x) {
                return(unlist(sapply(1:length(ind), function(grp) {
                  if (any(x == ind[[grp]])) return(grp) else return(NULL)
                })))
            })
            return(indBetas)
        }
        GetWeight <- function(p, indBetas, groupweights, groupnumber) {
            averageweight <- unlist(sapply(1:p, function(x) {
                grps <- indBetas[[x]]
                size <- length(grps)
                return(mean(groupweights[groupnumber \%in\% grps]))
            }))
            return(averageweight)
        }
        if (all(sapply(obtainHierarchy(Grouping), length) == 
            1)) 
            leaves <- 1:length(Grouping)
        else leaves <- which(sapply(obtainHierarchy(grouping.grouplvl), 
            length) == 1)
        IndBetas <- getIndBetas(p = length(values), Grouping)
        avrg <- sapply(Grouping[leaves], function(x) {
            mean(values[x])
        })
        medians <- sapply(Grouping[leaves], function(x) {
            median(values[x])
        })
        ord <- order(avrg)
        leaves <- leaves[ord]
        newLeaf <- as.factor(sapply(IndBetas, function(x) {
            which(leaves \%in\% x)
        }))
        originalLeaf <- as.factor(sapply(IndBetas, function(x) {
            leaves[which(leaves \%in\% x)]
        }))
        originalLeaf <- factor(originalLeaf, levels = leaves, 
            labels = leaves)
        dfBeta <- data.frame()
        for (l in unique(dfGrps$Fold)) {
            dfBeta2 <- data.frame(index = 1:length(values), Weight = GetWeight(length(values), 
                indBetas = IndBetas, groupweights = dfGrps$Group.weight[dfGrps$Fold == 
                  l], groupnumber = dfGrps$Group[dfGrps$Fold == 
                  l]), Continuous.values = values, newLeafGrp = newLeaf, 
                originalLeafGrp = originalLeaf)
            dfBeta2$Fold <- l
            dfBeta2$AverageGroupValue <- avrg[ord[dfBeta2$originalLeafGrp]]
            dfBeta2$MedianGroupValue <- medians[ord[dfBeta2$originalLeafGrp]]
            dfBeta <- rbind(dfBeta, dfBeta2)
        }
        dfBeta$Fold <- as.factor(dfBeta$Fold)
        SummdfBeta <- dfBeta \%>\% group_by(newLeafGrp) \%>\% summarise(meanValue = mean(Continuous.values), 
            minValue = min(Continuous.values), maxValue = max(Continuous.values), 
            q50Value = quantile(Continuous.values, 0.5), meanWeight = mean(Weight), 
            q25Weight = quantile(Weight, 0.25), q75Weight = quantile(Weight, 
                0.75), maxWeight = max(Weight), minWeight = min(Weight), 
            q50Weight = quantile(Weight, 0.5), originalLeafGrp = originalLeafGrp[1]) \%>\% 
            ungroup()
        if (any(is.na(SummdfBeta$meanValue))) {
            missingGroup <- which(is.na(SummdfBeta$meanValue))
            SummdfBeta[missingGroup, "minValue"] <- min(SummdfBeta$minValue, 
                na.rm = T)
            SummdfBeta[missingGroup, "maxValue"] <- max(SummdfBeta$maxValue, 
                na.rm = T)
            SummdfBeta[missingGroup, "meanValue"] <- mean(c(SummdfBeta[[missingGroup, 
                "minValue"]], SummdfBeta[[missingGroup, "maxValue"]]))
            levels(SummdfBeta$originalLeafGrp)[missingGroup] <- paste(levels(SummdfBeta$originalLeafGrp)[missingGroup], 
                " (missing data group)", sep = "")
        }
        if (!(boxplot | jitter)) {
            p1 <- ggplot(SummdfBeta) + geom_hline(yintercept = 1, 
                col = "grey", linetype = "dashed", size = ls, 
                alpha = 0.6) + geom_linerange(aes(x = q50Value, 
                ymax = q75Weight, ymin = q25Weight, col = originalLeafGrp), 
                linetype = "dashed") + geom_point(aes(x = q50Value, 
                y = minWeight, col = originalLeafGrp), shape = 2) + 
                geom_point(aes(x = q50Value, y = meanWeight, 
                  col = originalLeafGrp), shape = 1) + geom_point(aes(x = q50Value, 
                y = maxWeight, col = originalLeafGrp), shape = 6) + 
                geom_segment(aes(x = minValue, xend = maxValue, 
                  y = q50Weight, yend = q50Weight, col = originalLeafGrp)) + 
                scale_color_discrete(name = "Group") + labs(x = "Continuous co-data value", 
                y = "Prior variance weight") + theme_bw() + theme(axis.text.x = element_text(size = 12), 
                axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14), 
                axis.title.y = element_text(size = 14), legend.text = element_text(size = 12), 
                legend.title = element_text(size = 14))
        }
        else {
            dfBeta2 <- dfBeta \%>\% group_by(Fold) \%>\% distinct(originalLeafGrp, 
                .keep_all = T) \%>\% ungroup()
            p1 <- ggplot(dfBeta2) + geom_hline(yintercept = 1, 
                col = "grey", linetype = "dashed", size = ls, 
                alpha = 0.6) + geom_segment(data = SummdfBeta, 
                aes(x = minValue, xend = maxValue, y = q50Weight, 
                  yend = q50Weight, col = originalLeafGrp), size = ls) + 
                labs(x = "Continuous co-data value", y = "Prior variance weight") + 
                scale_color_discrete(name = "Group") + theme_bw() + 
                theme(axis.text.x = element_text(size = 12), 
                  axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14), 
                  axis.title.y = element_text(size = 14), legend.text = element_text(size = 12), 
                  legend.title = element_text(size = 14))
            if (boxplot) {
                p1 <- p1 + geom_boxplot(aes(x = MedianGroupValue, 
                  col = originalLeafGrp, y = Weight), width = widthBoxplot, 
                  varwidth = F, position = "identity", outlier.shape = NA)
            }
            if (jitter) {
                p1 <- p1 + geom_jitter(data = dfBeta2, aes(x = MedianGroupValue, 
                  col = originalLeafGrp, y = Weight), size = ps, 
                  alpha = 0.6, height = 0, width = widthBoxplot/4)
            }
        }
        return(p1)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
