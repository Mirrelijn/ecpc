\name{cv.ecpc}
\alias{cv.ecpc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Cross-validation for ecpc
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Cross-validates ecpc and returns model fit, summary statistics and cross-validated performance measures.
}
\usage{
cv.ecpc(Y,X,type.measure=c("MSE","AUC"),outerfolds=10,
        lambdas=NULL,ncores=1,balance=T,silent=F,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
Response data; n-dimensional vector (n: number of samples) for linear and logistic outcomes, or \code{\link{Surv}} object for Cox survival.
}
  \item{X}{
%%     ~~Describe \code{X} here~~
Observed data; (nxp)-dimensional matrix (p: number of covariates) with each row the observed high-dimensional feature vector of a sample.
}
  \item{type.measure}{
  Type of cross-validated performance measure returned.
}
  \item{outerfolds}{
  Number of cross-validation folds.
}
  \item{lambdas}{
  A vector of global ridge penalties for each fold; may be given, else estimated.
}
  \item{ncores}{
  Number of cores; if larger than 1, the outer cross-validation folds are processed in parallel over 'ncores' clusters.
}
  \item{balance}{
  (logistic, Cox) Should folds be balanced in response?
}
  \item{silent}{
%%     ~~Describe \code{normalise} here~~
Should output messages be suppressed (default FALSE)?
}
  \item{...}{
  Additional arguments used in \code{\link{ecpc}}.
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
A list with the following elements:
\item{ecpc.fit}{List with the ecpc model fit in each fold.}
\item{dfPred}{Data frame with information about out-of-bag predictions.}
\item{dfGrps}{Data frame with information about estimated group and grouping weights across folds.}
\item{dfCVM}{Data frame with cross-validated performance metric.}
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
#####################
# Simulate toy data #
#####################
p<-300 #number of covariates
n<-100 #sample size training data set
n2<-100 #sample size test data set

#simulate all betas i.i.d. from beta_k~N(mean=0,sd=sqrt(0.1)):
muBeta<-0 #prior mean
varBeta<-0.1 #prior variance
indT1<-rep(1,p) #vector with group numbers all 1 (all simulated from same normal distribution)

#simulate test and training data sets:
Dat<-simDat(n,p,n2,muBeta,varBeta,indT1,sigma=1,model='linear') 
str(Dat) #Dat contains centered observed data, response data and regression coefficients

##########################
# Make co-data groupings #
##########################
#Grouping 1: G random groups
G <- 5 #number of groups
#sample random categorical co-data:
categoricalRandom <- sample(1:G,p,replace=T) 
#make grouping, i.e. list with G groups:
groupingRandom <- lapply(1:G,function(g){which(categoricalRandom==g)}) 

#Grouping 2: informative hierarchical grouping
continuousCodata <- abs(Dat$beta) #use the magnitude of beta as continuous co-data
#Use adaptive discretisation to find a good discretisation of the continuous co-data;
# discretise in groups of covariates of various sizes:
groupingHierarchical <- splitMedian(values=continuousCodata,index = 1:p,
                        minGroupSize = 50,split="both") 
# and obtain grouping on group level that defines the hierarchy:
hierarchy.grouplevel <- obtainHierarchy(grouping = groupingHierarchical) 
#visualise hierarchical groups:
#visualiseGrouping(Grouping = groupingHierarchical,grouping.grouplvl = hierarchy.grouplevel) 
 
#######################
# Cross-validate ecpc #
#######################
# tic<-proc.time()[[3]]
# cv.fit <- cv.ecpc(type.measure="MSE",outerfolds=10,
#                   Y=Dat$Y,X=Dat$Xctd,
#                   groupings=list(groupingRandom,groupingHierarchical),
#                   groupings.grouplvl=list(NULL,hierarchy.grouplevel),
#                   hypershrinkage=c("ridge","hierLasso,ridge"),
#                   model="linear",maxsel=c(5,10,15,20))
# toc <- proc.time()[[3]]-tic
#
# str(cv.fit$ecpc.fit) #list containing the model fits on the folds
# str(cv.fit$dfPred) #data frame containing information on the predictions
# cv.fit$dfCVM #data frame with the cross-validated performance for ecpc with/without posterior selection and ordinary ridge

################
# Plot results #
################
# #Plot estimated grouping weights across folds
# visualiseGroupingweights(dfGrps = cv.fit$dfGrps,GroupingNames = c("random","hierarchical"),boxplot=F)
# 
# #Plot estimated group weights across folds
# #Grouping 1:
# visualiseGroupweights(dfGrps = cv.fit$dfGrps[cv.fit$dfGrps$Grouping==1,],Grouping=groupingRandom,boxplot = F)
# 
# #Grouping 2; plot group weight versus continuous co-data value:
# visualiseGroupweights(dfGrps = cv.fit$dfGrps[cv.fit$dfGrps$Grouping==2,],Grouping=groupingHierarchical,
#                       values=continuousCodata,grouping.grouplvl = hierarchy.grouplevel,boxplot = F)
# #visualise which hierarchical groups are selected in a specific fold:
# fold <- 1
# visualiseGrouping(Grouping = groupingHierarchical,grouping.grouplvl = hierarchy.grouplevel,
#                   groupweights = cv.fit$dfGrps$Group.weight[cv.fit$dfGrps$Grouping==2&cv.fit$dfGrps$Fold==fold]) 
}
