\name{visualiseGrouping}
\alias{visualiseGrouping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Visualise a grouping
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
visualiseGrouping(Grouping, groupweights, grouping.grouplvl, nodeSize = 10, ls = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Grouping}{
%%     ~~Describe \code{Grouping} here~~
}
  \item{groupweights}{
%%     ~~Describe \code{groupweights} here~~
}
  \item{grouping.grouplvl}{
%%     ~~Describe \code{grouping.grouplvl} here~~
}
  \item{nodeSize}{
%%     ~~Describe \code{nodeSize} here~~
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
function (Grouping, groupweights, grouping.grouplvl, nodeSize = 10, 
    ls = 1) 
{
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    cols <- c("#FC8D59", "ivory", "#91BFDB")
    Kg <- sapply(Grouping, length)
    G <- length(Grouping)
    groupNames <- names(Grouping)
    if (length(names(Grouping)) == 0) 
        groupNames <- 1:G
    if (missing(grouping.grouplvl) | length(grouping.grouplvl) == 
        0) {
        basicGraph <- T
        grouping.grouplvl <- lapply(1:G, function(x) x)
    }
    else basicGraph <- F
    i <- unlist(sapply(1:G, function(x) {
        rep(x, unlist(Kg)[x])
    }))
    j <- unlist(unlist(Grouping))
    ind <- sparseMatrix(i, j, x = 1)
    Zt <- ind
    if (G[1] > 1) {
        Zt[1:G[1], ] <- t(t(ind[1:G[1], ])/apply(ind[1:G[1], 
            ], 2, sum))
    }
    vertices <- data.frame(name = groupNames)
    if (missing(groupweights)) {
        vertices <- vertices \%>\% mutate(name2 = name, Selected = T, 
            GrpWeight = rep(1, G))
        showFillLegend <- F
    }
    else {
        vertices <- vertices \%>\% mutate(GrpWeight = round(groupweights, 
            digits = 2), Selected = groupweights != 0, name2 = paste(name, 
            " (", GrpWeight, ")", sep = ""))
        showFillLegend <- T
    }
    temp <- obtainHierarchy(grouping.grouplvl)
    children <- lapply(1:length(temp), function(i) {
        children <- setdiff(temp[[i]], i)
        for (j in children) {
            children <- setdiff(children, setdiff(temp[[j]], 
                j))
        }
        if (length(children) == 0) {
            return(NULL)
        }
        return(children)
    })
    edgeList <- data.frame(from = as.character(unlist(sapply(1:length(children), 
        function(x) {
            rep(x, length(children[[x]]))
        }))), to = as.character(unlist(children)))
    mygraph <- graph_from_data_frame(edgeList, vertices = vertices)
    if (basicGraph) {
        mygraph <- graph_from_data_frame(edgeList, vertices = vertices)
        lay <- create_layout(mygraph, layout = "linear")
        p <- ggraph(lay) + geom_node_label(aes(label = name2, 
            fill = GrpWeight, col = Selected), show.legend = showFillLegend, 
            fontface = "bold", size = nodeSize, repel = T, point.padding = NA, 
            box.padding = 0, force = 0.1) + scale_fill_gradientn(colors = c("white", 
            cols), limits = c(0, NA), values = scales::rescale(c(0, 
            min(10^-6, lay$GrpWeight[lay$GrpWeight > 0]), min(10^-6, 
                lay$GrpWeight[lay$GrpWeight > 0]), 1 - 1e-06, 
            1 - 1e-06, 1 + 1e-06, 1 + 1e-06, max(lay$GrpWeight))), 
            name = "Group weight") + scale_color_manual(values = c(`TRUE` = "black", 
            `FALSE` = cbPalette[1])) + theme_void()
        return(p)
    }
    nodesNoParents <- groupNames[!groupNames \%in\% edgeList$to]
    vertices <- rbind(vertices, rep(T, dim(vertices)[2]))
    vertices[dim(vertices)[1], c("name", "Selected")] <- c("origin", 
        T)
    vertices$visible <- T
    vertices$visible[vertices$name == "origin"] <- F
    edgeList <- rbind(edgeList, data.frame(from = rep("origin", 
        length(nodesNoParents)), to = as.character(nodesNoParents)))
    edgeList <- edgeList \%>\% mutate(visible = from != "origin")
    mygraph <- graph_from_data_frame(edgeList, vertices = vertices)
    lay <- create_layout(mygraph, layout = "dendrogram", circular = F)
    p <- ggraph(lay) + geom_edge_link(aes(alpha = node1.visible), 
        edge_width = ls, arrow = arrow(ends = "last", length = unit(3, 
            "mm"), type = "closed"), end_cap = circle(6, "mm")) + 
        scale_edge_alpha_discrete(range = c(0, 1), guide = F) + 
        geom_node_label(aes(label = name2, fill = GrpWeight, 
            col = Selected), show.legend = showFillLegend, fontface = "bold", 
            size = nodeSize) + scale_fill_gradientn(colors = c("white", 
        cols), limits = c(0, NA), values = scales::rescale(c(0, 
        min(lay$GrpWeight[lay$GrpWeight > 0], 10^-6), min(10^-6, 
            lay$GrpWeight[lay$GrpWeight > 0]), 1 - 1e-06, 1 - 
            1e-06, 1 + 1e-06, 1 + 1e-06, max(lay$GrpWeight))), 
        name = "Group weight") + scale_color_manual(values = c(`FALSE` = cbPalette[1], 
        `TRUE` = "black")) + coord_cartesian(ylim = c(min(lay$y) - 
        0.5, max(lay$y[lay$name != "origin"]) + 0.5), xlim = c(min(lay$x) - 
        0.5, max(lay$x) + 0.5)) + theme_void()
    return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
