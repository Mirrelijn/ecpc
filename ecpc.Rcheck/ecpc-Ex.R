pkgname <- "ecpc"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "ecpc-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('ecpc')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cv.ecpc")
### * cv.ecpc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cv.ecpc
### Title: Cross-validation for ecpc
### Aliases: cv.ecpc

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cv.ecpc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ecpc")
### * ecpc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ecpc
### Title: ecpc co-data learning
### Aliases: ecpc

### ** Examples

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
 
############ 
# Fit ecpc #
############
#fit ecpc for the two groupings, with ridge hypershrinkage for grouping 1, 
# and hierarchical lasso and ridge for grouping 2.
tic<-proc.time()[[3]]
fit <- ecpc(Y=Dat$Y,X=Dat$Xctd,groupings=list(groupingRandom,groupingHierarchical),
           groupings.grouplvl=list(NULL,hierarchy.grouplevel),
           hypershrinkage=c("ridge","hierLasso,ridge"),
           model="linear",maxsel=c(5,10,15,20),
           Y2=Dat$Y2,X2=Dat$X2ctd)
toc <- proc.time()[[3]]-tic

fit$tauglobal #estimated global prior variance
fit$gamma #estimated group weights (concatenated for the groupings)
fit$w #estimated grouping weights
summary(fit$beta) #estimated regression coefficients
summary(fit$betaPost) #estimated regression coefficients after posterior selection

c(fit$MSEecpc,fit$MSEridge) #mean squared error on test set for ecpc and ordinary ridge
fit$MSEPost #MSE on the test set of ecpc after posterior selection




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ecpc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("obtainHierarchy")
### * obtainHierarchy

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: obtainHierarchy
### Title: Obtain hierarchy
### Aliases: obtainHierarchy

### ** Examples

cont.codata <- seq(0,1,length.out=20) #continuous co-data
grouping <- splitMedian(values=cont.codata,split="lower",minGroupSize=5) #only split at lower continous co-data group
grouping.grouplvl <- obtainHierarchy(grouping) #obtain groups on group level defining the hierarchy




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("obtainHierarchy", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simDat")
### * simDat

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simDat
### Title: Simulate data
### Aliases: simDat

### ** Examples

# #required libraries
# library(mvtnorm)
#
# n<-10
# p<-30
# #simulate beta from two normal distributions; beta_k ~ N(mu_k,tau^2_k)
# muGrp <- c(0,0.1) #mean (mu_1,mu_2)
# varGrp <- c(0.05,0.01) #variance (tau^2_1,tau^2_2)
# indT <- rep(c(1,2),each=15) #group number of each covariate; first half in group 1, second half in group 2
#
# dataLin <- simDat(n, p, n2 = 20, muGrp, varGrp, indT, sigma = 1, model = "linear", 
#     flag = T)
# dataLog <- simDat(n, p, n2 = 20, muGrp, varGrp, indT, model = "logistic", 
#     flag = T)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simDat", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("splitMedian")
### * splitMedian

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: splitMedian
### Title: Discretise continuous data in multiple granularities
### Aliases: splitMedian

### ** Examples


cont.codata <- seq(0,1,length.out=20) #continuous co-data
grouping1 <- splitMedian(values=cont.codata,minGroupSize=5) #full tree with minimum group size 5
grouping2 <- splitMedian(values=cont.codata,split="lower",minGroupSize=5) #only split at lower continous co-data group

part <- sample(1:length(cont.codata),15) #discretise only for a part of the continuous co-data
cont.codata[-part] <- NaN #suppose rest is missing
grouping3 <- splitMedian(values=cont.codata[part],index=part,minGroupSize=5) #make grouping of non-missing values
grouping3 <- c(grouping3,list(which(is.nan(cont.codata)))) #add missing data group





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("splitMedian", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("visualiseGrouping")
### * visualiseGrouping

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: visualiseGrouping
### Title: Visualise a grouping
### Aliases: visualiseGrouping

### ** Examples

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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("visualiseGrouping", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("visualiseGroupingweights")
### * visualiseGroupingweights

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: visualiseGroupingweights
### Title: Visualise estimated grouping weights
### Aliases: visualiseGroupingweights

### ** Examples

# #required libraries
# library(ggplot2)
# library(dplyr)
#
# dfGrps <- data.frame(Grouping=rep(c(1,2),each=10),
#                      Grouping.weight=c(rnorm(10,0,0.01),rnorm(10,1,0.05)),
#                      Fold=rep(1:10,2))
# GroupingNames <- c("Grouping1","Grouping2")
# visualiseGroupingweights(dfGrps, GroupingNames, hist = F, boxplot = T,jitter=T)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("visualiseGroupingweights", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("visualiseGroupweights")
### * visualiseGroupweights

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: visualiseGroupweights
### Title: Visualise estimated group weights
### Aliases: visualiseGroupweights

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("visualiseGroupweights", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
