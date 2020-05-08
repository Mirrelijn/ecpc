#Demo ecpc
#NOTE: in Rstudio; press ALT+o for overview of code section headings
#load libraries
library(mvtnorm) #for simulating data set with rmvnorm
library(ggplot2) #for plots
library(igraph) #for plots
library(ggraph) #for plots
library(dplyr) #for data handling
library(gglasso) #for hierarchical groups
library(Matrix)
library(glmnet)
library(penalized)
library(foreach) #for parallel computations
library(doParallel) #for parallel computations

#Download package from github;
#library(devtools); install_github("Mirrelijn/ecpc/Rpackage")
library(ecpc)

#Simulate data----
p<-300 #number of covariates
n<-100 #sample size training data set
n2<-100 #sample size test data set
#simulate betas i.i.d. from beta_k~N(mean=0,sd=sqrt(0.1))
muBeta<-0 #prior mean
varBeta<-0.1 #prior variance
indT1<-rep(1,p) #vector with group numbers; simulate each beta from the same normal distribution
Dat<-simDat(n,p,n2,muBeta,varBeta,indT1,sigma=1,model='linear') #simulate test and training data sets
str(Dat) #Dat contains centered observed data, response data and the true vector of regression coefficients

#Make co-data groupings----
G <- 5 #number of random/informative groups

#Grouping 1: G random non-overlapping groups
categoricalRandom <- sample(1:G,p,replace=T) #simulate random categorical co-data
groupingRandom <- lapply(1:G,function(g){which(categoricalRandom==g)}) #make grouping; list with G groups

#Grouping 2: G informative non-overlapping groups of the same size as the random groups
rankBeta<-order(abs(Dat$beta)) #betas ranked in order of magnitude
cumgroupSizes <- c(1,cumsum(sapply(groupingRandom,length))) #cumulated group sizes
groupingInformative <- lapply(1:G,function(g){rankBeta[cumgroupSizes[g]:cumgroupSizes[g+1]]}) #make informative grouping
#sapply(groupingInformative,function(group){mean(abs(Dat$beta[group]))}) #check if magnitude indeed increases

#Grouping 3: informative hierarchical grouping based on magnitude beta
continuousCodata <- abs(Dat$beta) #use the magnitude of beta as continuous co-data
#Use adaptive discretisation to find a good discretisation of the continuous co-data;
groupingHierarchical <- splitMedian(values=continuousCodata,index = 1:p,minGroupSize = 50,split="both") #discretise in groups of covariates
hierarchy.grouplevel <- obtainHierarchy(grouping = groupingHierarchical) #obtain grouping on group level that defines the hierarchy
visualiseGrouping(Grouping = groupingHierarchical,grouping.grouplvl = hierarchy.grouplevel) #visualise hierarchical groups

#Fit ecpc----
#fit ecpc for the three groupings, with ridge hypershrinkage for grouping 1 and 2, and hierarchical Lasso with ridge for the
#hierarchical groups corresponding to the adaptive discretisation of the continuous co-data.
#Takes ~1 minute to fit
tic<-proc.time()[[3]]
fit <- ecpc(Y=Dat$Y,X=Dat$Xctd,groupings=list(groupingRandom,groupingInformative,groupingHierarchical),
            groupings.grouplvl=list(NULL,NULL,hierarchy.grouplevel),
            hypershrinkage=c("ridge","ridge","hierLasso,ridge"),
            model="linear",maxsel=c(5,10,15,20),
            Y2=Dat$Y2,X2=Dat$X2ctd)
toc <- proc.time()[[3]]-tic

summary(fit$beta) #beta contains the estimated regression coefficients
summary(fit$betaPost) #betaPost contains the estimated regression coefficients after posterior selection
c(fit$MSEecpc,fit$MSEridge) #MSE on the test set
fit$MSEPost #MSE of the parsimonious models on the test set

fit$w #grouping weights
fit$gamma #group weights are concatenated for the groupings
fit$tauglobal #global prior variance

#Cross-validate ecpc----
#Takes ~'outerfolds' minutes
folds2 <- produceFolds(n,3,Dat$Y,model="linear")
tic<-proc.time()[[3]]
cv.fit <- cv.ecpc(type.measure="MSE",outerfolds=folds2,
                  ncores=1, Y=Dat$Y,X=Dat$Xctd,
                  groupings=list(groupingRandom,groupingInformative,groupingHierarchical),
                  groupings.grouplvl=list(NULL,NULL,hierarchy.grouplevel),
                  hypershrinkage=c("ridge","ridge","hierLasso,ridge"),
                  model="linear",maxsel=c(5,10,15,20))
toc <- proc.time()[[3]]-tic

# #or execute in parallel, use detectCores() to view number of available cores
#Takes ~'outerfolds'/'ncores' 
tic<-proc.time()[[3]]
cv.fit<- cv.ecpc(type.measure="MSE",outerfolds=10,
                 ncores=detectCores(), Y=Dat$Y,X=Dat$Xctd,
                 groupings=list(groupingRandom,groupingInformative,groupingHierarchical),
                 groupings.grouplvl=list(NULL,NULL,hierarchy.grouplevel),
                 hypershrinkage=c("ridge","ridge","hierLasso,ridge"),
                 model="linear",maxsel=c(5,10,15,20))
toc <- proc.time()[[3]]-tic

str(cv.fit$ecpc.fit) #list containing the model fits on the folds
str(cv.fit$dfPred) #data frame containing information about the predictions
cv.fit$dfCVM #data frame with the cross-validated performance for ecpc with/without posterior selection and ordinary ridge
cv.fit$dfCVM %>% group_by(Method) %>% summarise(mean(CVM)) #for MSE; compute mean MSE across folds

#Plot results----
#Plot estimated grouping weights across folds
visualiseGroupingweights(dfGrps = cv.fit$dfGrps,GroupingNames = c("random","informative","hierarchical"),boxplot=F)

#Plot estimated group weights across folds
#Grouping 1:
visualiseGroupweights(dfGrps = cv.fit$dfGrps[cv.fit$dfGrps$Grouping==1,],Grouping=groupingRandom,boxplot = F)

#Grouping 2:
visualiseGroupweights(dfGrps = cv.fit$dfGrps[cv.fit$dfGrps$Grouping==2,],Grouping=groupingInformative,boxplot = F)

#Grouping 3; plot group weight versus continuous co-data value:
visualiseGroupweights(dfGrps = cv.fit$dfGrps[cv.fit$dfGrps$Grouping==3,],Grouping=groupingHierarchical,
                      values=continuousCodata,grouping.grouplvl = hierarchy.grouplevel,boxplot = F)
#visualise which hierarchical groups are selected in a specific fold:
fold <- 1
visualiseGrouping(Grouping = groupingHierarchical,grouping.grouplvl = hierarchy.grouplevel,
                  groupweights = cv.fit$dfGrps$Group.weight[cv.fit$dfGrps$Grouping==3&cv.fit$dfGrps$Fold==fold])



