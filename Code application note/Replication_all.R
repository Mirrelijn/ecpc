#Replication script for all analyses, figures and tables in the application note for ecpc
#Note: press alt+O in Rstudio to fold code to headings

#setwd("...") #set working directory to desired folder
#knitr::spin("Replication_all.R")
#make sure the following packages are installed correctly:
if(!requireNamespace("ggplot2")) install.packages("ggplot2")
if(!requireNamespace("ggpubr")) install.packages("ggpubr")
if(!requireNamespace("ecpc")) install.packages("ecpc")
if(!requireNamespace("squeezy")) install.packages("squeezy")
if(!requireNamespace("dplyr")) install.packages("dplyr")
if(!requireNamespace("ggplot2")) install.packages("ggplot2")
if(!requireNamespace("RColorBrewer")) install.packages("RColorBrewer")
if(!requireNamespace("foreach")) install.packages("foreach")
if(!requireNamespace("doParallel")) install.packages("doParallel")
if(!requireNamespace("glmnet")) install.packages("glmnet")
if(!requireNamespace("mvtnorm")) install.packages("mvtnorm")
if(!requireNamespace("pROC")) install.packages("pROC")
if(!requireNamespace("devtools")) install.packages("devtools")
library(devtools)
if(!requireNamespace("fwelnet")) install_github("kjytay/fwelnet")
if(!requireNamespace("CoRF")) install_github("DennisBeest/CoRF") 


#Default settings to run everything----
#Section 3.1 settings
runLinear <- TRUE
runGAM <- TRUE
runSCAM <- TRUE
runSCAMpmi <- TRUE
runAD <- TRUE

generateData <- TRUE #set to true if not generated before
nSim <- 50 #set to lower number for quicker run
runParallel <- TRUE

#section 3.2 settings 
run_ecpc <- TRUE
run_glmnet <- TRUE
run_fwelnet <- TRUE

#Section 4
run_example <- TRUE #example code in the paper

#full analyses for 4 different settings:
run_withConstraints <- TRUE
run_withoutConstraints <- TRUE
run_transformed <- TRUE
run_withConstraints2 <- TRUE




###############################################################################

#####Section 1&2: short examples###############################################
###############################################################################
figno <- 1
wdthpdf <- 8
hghtpdf <- 5

library(ggplot2) #for saving plots
library(ggpubr)

#Section 1.1----
#load package
library(ecpc)

#section 2.1----
#Simulate linear response data----
set.seed(1)
p<-300 #number of covariates
n<-100 #sample size training data set
n2<-100 #sample size test data set
beta <- rnorm(p, mean=0, sd=0.1) #simulate effects
X <- matrix(rnorm(n*p,mean=0,sd=1), n,p) #simulate observed training data
Y <- rnorm(n,mean = X%*%beta, sd = 1) #simulate response training data
X2 <- matrix(rnorm(n2*p,mean=0,sd=1), n,p) #simulate observed test data
Y2 <- rnorm(n2,mean = X2%*%beta, sd = 1) #simulate response test data


#Consider two co-data variables----
Z1 <- abs(beta) #informative co-data
Z2 <- rnorm(p,mean=0,sd=1) #random, non-informative co-data
Z <- cbind( Z1 , Z2 ) #(px2)-dimensional co-data matrix

#Fit model----
fit <- ecpc(Y,X,Z=list(Z),X2=X2,Y2=Y2)
summary(fit$beta) #fitted regression coefficients
fit$gamma*fit$tauglobal #fitted co-data weights

fit$MSEecpc
fit$MSEridge

print(fit) #print the estimated prior variables
summary(fit) #print summaries of the estimated regression coefficients and prior variances

p1 <- plot(fit, show="coefficients") #plot regression coefficients
p1
p2 <- plot(fit, show="priorweights", Z=list(Z)) #plot prior variances per co-data variable
p2
ggarrange(plotlist=list(p1,p2))
ggsave(filename=paste("FigExample",figno,".pdf",sep=""),
       width=1.5*wdthpdf, height=hghtpdf)
figno <- figno + 1

new_penalties <- penalties(fit, tauglobal = fit$tauglobal * 2, Z=list(Z)) #compute penalties     
# corresponding to prior variances twice as large
new_coefficients <- coef(fit, penalties=new_penalties, X=X, Y=Y) #estimate corresponding
# regression coefficients

#compute penalties for given prior parameters 
new_penalties2 <- penalties(tauglobal = 1, sigmahat = 1, gamma = c(1,1), w = 1, Z=list(Z)) 
#skip prior parameter estimation and estimate regression coefficients only
new_coefficients2 <- coef.ecpc(penalties=new_penalties2, X=X, Y=Y) 

#Section 2.2----
#Create spline basis matrix for the two co-data variables----
Z1.s <- createZforSplines(values=Z1, G=20, bdeg=3) #use 20 B-splines of 3rd degree
S1.Z1 <- createS(orderPen=2, G=20) #create 2nd order difference penalty matrix
Z2.s <- createZforSplines(values=Z2, G=30, bdeg=3) #use 30 B-splines functions
S1.Z2 <- createS(orderPen=2, G=30)

Z.all <- list(Z1=Z1.s,Z2=Z2.s)
paraPen.all <- list(Z1=list(S1=S1.Z1), Z2=list(S1=S1.Z2))

#Fit model with splines----
fit.gam <- ecpc(Y, X, Z = Z.all, paraPen = paraPen.all, intrcpt.bam=TRUE, X2=X2, Y2=Y2)
head(fit.gam$gamma*fit.gam$tauglobal) #fitted co-data spline coefficients
fit.gam$gamma0*fit.gam$tauglobal #fitted co-data intercept 

p1<-plot(fit.gam, show="priorweights", Z=Z.all) #co-data variables on the x-axis
p1
values <- list(Z1, Z2)
p2<-plot(fit.gam, show="priorweights", Z=Z.all, values = values) #continuous values the x-axis
p2
ggarrange(plotlist=list(p1,p2))
ggsave(filename=paste("FigExample",figno,".pdf",sep=""),
       width=3*wdthpdf, height=hghtpdf)
figno <- figno + 1

#plot contribution of one co-data source directly----
i <- 2 #1 for informative, 2 for non-informative
codataNO <- attributes(fit.gam$gamma)$codataSource #corresponding co-data source number in fit.gam$gamma
sk <- as.vector(Z.all[[i]]%*%fit.gam$gamma[codataNO==i])*fit.gam$tauglobal
par(mfrow=c(1,1))
plot(Z[,i],sk)

#Section 2.3----
#Create constraints for the two co-data variables----
Con.Z1 <- createCon(G=20, shape="positive+monotone.i") 
Con.Z2 <- createCon(G=30, shape="convex") 
paraCon <- list(Z1=Con.Z1, Z2=Con.Z2)

#Fit model with shape constrained splines----
fit.scam <- ecpc(Y, X, Z = Z.all, paraPen = paraPen.all, paraCon = paraCon, X2=X2, Y2=Y2)
fit.scam$gamma*fit.scam$tauglobal #fitted co-data spline coefficients
fit.scam$w #fitted co-data weights
fit.scam$gamma0*fit.scam$tauglobal #fitted co-data intercept excluded by default

p1<-plot(fit.scam, show="priorweights", Z=Z.all, values=values)
ggarrange(plotlist=list(p1))
ggsave(filename=paste("FigExample",figno,".pdf",sep=""),
       width=1.5*wdthpdf, height=hghtpdf)
figno <- figno + 1

#plot shape-constrained contribution
i <- 2 #1 for informative, 2 for non-informative
codataNO <- c(rep(1,20), rep(2,30)) #corresponding groupset number in fit$gamma
sk <- as.vector(Z.all[[i]]%*%fit.scam$gamma[codataNO==i])*fit.scam$tauglobal
par(mfrow=c(1,1))
plot(Z[,i],sk)

#Section 2.4----
#Transform to elastic net penalties----
library(squeezy)

#Use squeezy function to transform estimated ridge penalties to elastic net penalties
#for some alpha, here for alpha=0.3
fit.EN <- squeezy(Y,X,alpha=0.3, X2=X2, Y2=Y2,lambdas = fit$penalties)
summary(fit.EN$lambdapApprox) #transformed elastic net penalties
summary(fit.EN$betaApprox) #fitted elastic net regression coefficients






###############################################################################

#Section 3: simulation study####################################################
################################################################################
#Section 3.1 Simulation study comparing different co-data models----

#Load libraries-----
library(ecpc)
library(dplyr) #for data wrangling results
library(ggplot2) #for plotting results
library(RColorBrewer) #for plotting results
library(foreach) #for parallel computing
library(doParallel) #for parallel computing


if(runParallel){ #set to 1 to setup parallel backend to use many processors
  cores=detectCores()
  if(!("cl"%in%ls())){
    cl <- makeCluster(cores-1) #not to overload your computer
    registerDoParallel(cl)
  }
}


# Simulate toy data ------
p<-300 #number of covariates
n<-100 #sample size training data set
n2<-100 #sample size test data set
sigma <- 1
tauglobal <- 0.1 #prior variance

#simulate all betas i.i.d. from beta_k~N(mean=0,sd=sqrt(0.1)):
set.seed(7474)
Dat <- list()
Dat$beta <- rnorm(p,sd=sqrt(tauglobal))
Dat$Xctd <- matrix(rnorm(n*p) ,n,p)
Dat$X2ctd <- matrix(rnorm(n*p) ,n,p)
Dat$Y <- rnorm(n,mean=c(Dat$Xctd%*%Dat$beta),sd=sigma)

if(generateData){
  AllData <- list()
  for(i in 1:nSim){
    AllData[[i]] <- list()
    AllData[[i]]$beta <- Dat$beta
    AllData[[i]]$Xctd <- matrix(rnorm(n*p) ,n,p)
    AllData[[i]]$X2ctd <- matrix(rnorm(n*p) ,n,p)
    AllData[[i]]$Y <- rnorm(n,mean=c(AllData[[i]]$Xctd%*%AllData[[i]]$beta),sd=sigma)
    AllData[[i]]$Y2 <- rnorm(n2,mean=c(AllData[[i]]$X2ctd%*%AllData[[i]]$beta),sd=sigma)
  }
  save(Dat,AllData,file="SimData.Rdata")
}
load("SimData.Rdata")

#Co-data settings----
G <- c(20,50) #number of splines
Z.all <- list()
ZI.all <- list() #co-data matrix with intercept
Zs.all <- list()
S1 <- list()
Con.p <- list() #positivity constraints
Con.pmi <- list() #positivity and monotonically increasing constraints
#setting 1: non-informative, unequally spaced
Z.all[["noninformative"]] <- rnorm(p) #for linear co-data model
ZI.all[["noninformative"]] <- cbind(rep(1,p), Z.all[["noninformative"]])
#for generalised additive co-data model:
for(g in 1:length(G)){
  Zs.all[[g]] <- list()
  Zs.all[[g]][["noninformative"]] <- createZforSplines(values=Z.all[["noninformative"]], G=G[g], bdeg=3) 
  S1[[g]] <- createS(orderPen=2, G=G[g]) #create 2nd order difference penalty matrix (same for all co-data)
  
  Con.p[[g]] <- createCon(G=G[g], shape="positive") #create constraints
  Con.pmi[[g]] <- createCon(G=G[g], shape="positive+monotone.i") #create constraints
}
plot(Z.all[[1]],Dat$beta^2)

#setting 2: informative, unequally spaced, information at edge
Z.all[["size.edge"]] <- abs(Dat$beta)
ZI.all[["size.edge"]] <- cbind(rep(1,p), Z.all[["size.edge"]])
#for generalised additive co-data model:
for(g in 1:length(G)){
  Zs.all[[g]][["size.edge"]] <- createZforSplines(values=Z.all[["size.edge"]], G=G[g], bdeg=3) 
}
plot(Z.all[[2]],Dat$beta^2)

#setting 3: non-informative, transformed for equally spacing
Z.all[["noninformative+transformed"]] <- order(order(Z.all[["noninformative"]],
                                                     decreasing=FALSE),decreasing=FALSE)
ZI.all[["noninformative+transformed"]] <- cbind(rep(1,p), Z.all[["noninformative+transformed"]])
#for generalised additive co-data model:
for(g in 1:length(G)){
  Zs.all[[g]][["noninformative+transformed"]] <- createZforSplines(values=Z.all[["noninformative+transformed"]], G=G[g], bdeg=3) 
}
plot(Z.all[[3]],Dat$beta^2)

#setting 4: informative, transformed for equally spacing
Z.all[["size.edge+transformed"]] <- order(order(abs(Dat$beta),decreasing=FALSE),decreasing=FALSE)
ZI.all[["size.edge+transformed"]] <- cbind(rep(1,p), Z.all[["size.edge+transformed"]])
for(g in 1:length(G)){
  Zs.all[[g]][["size.edge+transformed"]] <- createZforSplines(values=Z.all[["size.edge+transformed"]], G=G[g], bdeg=3) 
}
plot(Z.all[[4]],Dat$beta^2)

#Compare with adaptive discretisation of continuous co-data----
groupsets <- list()
hierarchy.grouplevel <- list()
Z.AD <- list() #useful for computing prior variances

#setting 1: non-informative hierarchical discretisation
#Use adaptive discretisation to find a good discretisation of the continuous co-data;
# discretise in groups of covariates of various sizes:
groupsets[["noninformative"]] <- splitMedian(values=Z.all[["noninformative"]],index = 1:p,
                                             minGroupSize = 50,split="both") 
Z.AD[["noninformative"]] <- createZforGroupset(groupsets[["noninformative"]])
# and obtain group set on group level that defines the hierarchy:
hierarchy.grouplevel[["noninformative"]] <- obtainHierarchy(groupset = groupsets[["noninformative"]])

#setting 2: informative hierarchical discretisation
#Use adaptive discretisation to find a good discretisation of the continuous co-data;
# discretise in groups of covariates of various sizes:
groupsets[["size.edge"]] <- splitMedian(values=Z.all[["size.edge"]],index = 1:p,
                                        minGroupSize = 50,split="both") 
Z.AD[["size.edge"]] <- createZforGroupset(groupsets[["size.edge"]])
# and obtain group set on group level that defines the hierarchy:
hierarchy.grouplevel[["size.edge"]] <- obtainHierarchy(groupset = groupsets[["noninformative"]])


#Simulation settings----
logname <- "logSim.txt"
write(paste(Sys.time(),logname),file=logname,append=T)

#Fit ecpc on simulated data sets: linear co-data model----
fname <- paste(pathResults,"SimResAppNoteLinear",".Rdata",sep="") #(with intercept in ZI.all)

if(runLinear){
  #df <- data.frame()
  #for(sim in 1:nSim){
  finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
                         .packages = c("ecpc")) %dopar% {
                           df <- data.frame()
                           for(setting in 1:4){
                             for(method in "none"){#c("ML","fREML","GCV.Cp")){
                               tic<-proc.time()[[3]]
                               fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
                                           Z=ZI.all[setting],
                                           #bam.method=method,intrcpt.bam = F,
                                           model="linear",maxsel=c(5,10,15,20),
                                           Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,est_beta_method = "multiridge")
                               toc <- proc.time()[[3]]-tic
                               
                               vk <- (ZI.all[[setting]]%*%fit$gamma )*fit$tauglobal
                               #vk <- (matrix(Z.all[[setting]],p,1)%*%fit$gamma )*fit$tauglobal #without intercept
                               #plot(Z.all[[setting]],vk)
                               
                               
                               temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
                               temp$Z <- Z.all[[setting]]
                               temp$Time <- toc
                               temp$Covariate <- 1:p
                               temp$G <- 2
                               temp$setting <- setting
                               temp$bam.method <- method
                               temp$method <- "linear"
                               temp$MSEridge <- fit$MSEridge
                               temp$MSEecpc <- fit$MSEecpc
                               temp$Sim <- sim
                               temp$Codata <- names(Z.all)[setting]
                               temp$Transform <- F; if(setting>=3) temp$Transform <- T
                               temp$ZUntransformed <- Z.all[[setting]]
                               if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
                               df <- rbind(df,temp)
                             }
                           }
                           write(paste(Sys.time(),"simulation",sim,"of",nSim,"done"),file=logname,append=T)
                           list("df"=df)
                         }                                     
  #str(finalMatrix)
  df <- finalMatrix[[1]]; for(i in 2:nSim) df <- rbind(df,finalMatrix[[i]])
  save(finalMatrix,df,file=fname)
}


#Fit ecpc on simulated data sets: generalised additive co-data model----
fname <- paste(pathResults,"SimResAppNoteGAM.Rdata",sep="")
print(fname)

if(runGAM){
  #df <- data.frame()
  #for(sim in 1:nSim){
  finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
                         .packages = c("ecpc")) %dopar% {
                           df <- data.frame()
                           for(setting in 1:4){
                             for(g in 1:length(G)){
                               for(method in c("splits","ML","fREML","GCV.Cp")){
                                 tic<-proc.time()[[3]]
                                 
                                 if(method=="splits"){
                                   fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
                                               Z=Zs.all[[g]][setting],
                                               paraPen=list(Z1=list(S1=S1[[g]])),
                                               intrcpt.bam = F,
                                               model="linear",maxsel=c(5,10,15,20),
                                               Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
                                               hypershrinkage = "ridge",
                                               est_beta_method = "multiridge")
                                 }else{
                                   fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
                                               Z=Zs.all[[g]][setting],
                                               paraPen=list(Z1=list(S1=S1[[g]])),
                                               bam.method=method,intrcpt.bam = F,
                                               model="linear",maxsel=c(5,10,15,20),
                                               Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
                                               est_beta_method = "multiridge")
                                 }
                                 
                                 toc <- proc.time()[[3]]-tic
                                 
                                 vk <- (Zs.all[[g]][[setting]]%*%fit$gamma )*fit$tauglobal
                                 #plot(Z.all[[setting]],vk)
                                 
                                 
                                 
                                 temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
                                 temp$Z <- Z.all[[setting]]
                                 temp$Time <- toc
                                 temp$Covariate <- 1:p
                                 temp$G <- G[g]
                                 temp$setting <- setting
                                 temp$bam.method <- method
                                 temp$method <- "gam"
                                 temp$MSEridge <- fit$MSEridge
                                 temp$MSEecpc <- fit$MSEecpc
                                 temp$Sim <- sim
                                 temp$Codata <- names(Z.all)[setting]
                                 temp$Transform <- F; if(setting>=3) temp$Transform <- T
                                 temp$ZUntransformed <- Z.all[[setting]]
                                 if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
                                 df <- rbind(df,temp)
                               }
                             }
                           }
                           write(paste(Sys.time(),"simulation",sim,"of",nSim,"done"),file=logname,append=T)
                           list("df"=df)
                         }                                     
  #str(finalMatrix)
  df <- finalMatrix[[1]]; for(i in 2:nSim) df <- rbind(df,finalMatrix[[i]])
  save(finalMatrix,df,file=fname)
}

#Fit ecpc on simulated data sets: positive constrained generalised additive co-data model----
fname <- paste(pathResults,"SimResAppNoteSCAMp.Rdata",sep="")
print(fname)

if(runSCAM){
  #df <- data.frame()
  #for(sim in 1:nSim){
  finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
                         .packages = c("ecpc")) %dopar% {
                           df <- data.frame()
                           for(setting in 1:4){
                             for(g in 1:length(G)){
                               for(method in c("splits")){
                                 tic<-proc.time()[[3]]
                                 
                                 fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
                                             Z=Zs.all[[g]][setting],
                                             paraPen=list(Z1=list(S1=S1[[g]])),
                                             paraCon = list(Z1=Con.p[[g]]),
                                             intrcpt.bam = F,
                                             model="linear",maxsel=c(5,10,15,20),
                                             Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
                                             est_beta_method = "multiridge")
                                 toc <- proc.time()[[3]]-tic
                                 
                                 vk <- (Zs.all[[g]][[setting]]%*%fit$gamma )*fit$tauglobal
                                 plot(Z.all[[setting]],vk)
                                 
                                 
                                 
                                 temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
                                 temp$Z <- Z.all[[setting]]
                                 temp$Time <- toc
                                 temp$Covariate <- 1:p
                                 temp$G <- G[g]
                                 temp$setting <- setting
                                 temp$bam.method <- method
                                 temp$method <- "scam.p"
                                 temp$MSEridge <- fit$MSEridge
                                 temp$MSEecpc <- fit$MSEecpc
                                 temp$Sim <- sim
                                 temp$Codata <- names(Z.all)[setting]
                                 temp$Transform <- F; if(setting>=3) temp$Transform <- T
                                 temp$ZUntransformed <- Z.all[[setting]]
                                 if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
                                 df <- rbind(df,temp)
                               }
                             }
                           }
                           write(paste(Sys.time(),"simulation",sim,"of",nSim,"done"),file=logname,append=T)
                           list("df"=df)
                         }                                     
  #str(finalMatrix)
  df <- finalMatrix[[1]]; for(i in 2:nSim) df <- rbind(df,finalMatrix[[i]])
  save(finalMatrix,df,file=fname)
}




#Fit ecpc on simulated data sets: positive+monotone increasing constrained generalised additive co-data model----
fname <- paste(pathResults,"SimResAppNoteSCAMpmi.Rdata",sep="")
print(fname)

if(runSCAMpmi){
  #df <- data.frame()
  #for(sim in 1:nSim){
  finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
                         .packages = c("ecpc")) %dopar% {
                           df <- data.frame()
                           for(setting in 1:4){
                             for(g in 1:length(G)){
                               for(method in c("splits")){
                                 tic<-proc.time()[[3]]
                                 
                                 fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
                                             Z=Zs.all[[g]][setting],
                                             paraPen=list(Z1=list(S1=S1[[g]])),
                                             paraCon = list(Z1=Con.pmi[[g]]),
                                             intrcpt.bam = F,
                                             model="linear",maxsel=c(5,10,15,20),
                                             Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
                                             est_beta_method = "multiridge")
                                 toc <- proc.time()[[3]]-tic
                                 
                                 vk <- (Zs.all[[g]][[setting]]%*%fit$gamma )*fit$tauglobal
                                 #plot(Z.all[[setting]],vk)
                                 
                                 
                                 
                                 temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
                                 temp$Z <- Z.all[[setting]]
                                 temp$Time <- toc
                                 temp$Covariate <- 1:p
                                 temp$G <- G[g]
                                 temp$setting <- setting
                                 temp$bam.method <- method
                                 temp$method <- "scam.pmi"
                                 temp$MSEridge <- fit$MSEridge
                                 temp$MSEecpc <- fit$MSEecpc
                                 temp$Sim <- sim
                                 temp$Codata <- names(Z.all)[setting]
                                 temp$Transform <- F; if(setting>=3) temp$Transform <- T
                                 temp$ZUntransformed <- Z.all[[setting]]
                                 if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
                                 df <- rbind(df,temp)
                               }
                             }
                           }
                           write(paste(Sys.time(),"simulation",sim,"of",nSim,"done"),file=logname,append=T)
                           list("df"=df)
                         }                                     
  #str(finalMatrix)
  df <- finalMatrix[[1]]; for(i in 2:nSim) df <- rbind(df,finalMatrix[[i]])
  save(finalMatrix,df,file=fname)
}



#Fit ecpc on simulated data sets: adaptive discretisation co-data model----
fname <- paste(pathResults,"SimResAppNoteAD.Rdata",sep="")
print(fname)

if(runAD){
  #df <- data.frame()
  #for(sim in 1:nSim){
  finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
                         .packages = c("ecpc")) %dopar% {
                           df <- data.frame()
                           for(setting in 1:2){
                             #for(g in 1:length(G)){
                             for(method in c("splits")){
                               tic<-proc.time()[[3]]
                               
                               fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
                                           groupsets=groupsets[setting],
                                           groupsets.grouplvl = hierarchy.grouplevel[setting],
                                           hypershrinkage = "hierLasso,ridge",
                                           model="linear",maxsel=c(5,10,15,20),
                                           Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
                                           est_beta_method = "multiridge")
                               toc <- proc.time()[[3]]-tic
                               
                               vk <- as.vector(Z.AD[[setting]]%*%fit$gamma*fit$tauglobal)
                               #plot(Z.all[[setting]],vk)
                               
                               
                               temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
                               temp$Z <- Z.all[[setting]]
                               temp$Time <- toc
                               temp$Covariate <- 1:p
                               temp$G <- length(groupsets[[setting]])
                               temp$setting <- setting
                               temp$bam.method <- method
                               temp$method <- "AD"
                               temp$MSEridge <- fit$MSEridge
                               temp$MSEecpc <- fit$MSEecpc
                               temp$Sim <- sim
                               temp$Codata <- names(Z.all)[setting]
                               temp$Transform <- F; if(setting>=3) temp$Transform <- T
                               temp$ZUntransformed <- Z.all[[setting]]
                               if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
                               df <- rbind(df,temp)
                             }
                             #}
                           }
                           write(paste(Sys.time(),"simulation",sim,"of",nSim,"done"),file=logname,append=T)
                           list("df"=df)
                         }                                     
  #str(finalMatrix)
  df <- finalMatrix[[1]]; for(i in 2:nSim) df <- rbind(df,finalMatrix[[i]])
  save(finalMatrix,df,file=fname)
}




#Plots: general parameters----
wdth<-600
hght<-wdth*5/8
wdthpdf <- wdth/75
hghtpdf <- hght/75
ts <- 16 #basis text size in figures
ls <- 1.5 #basis line size in figures
ps <- 2 #basis point size in figures
sz <- 2 #point size
strk <- 1.5 #stroke size
palette <- "Dark2"
display.brewer.all(10,colorblindFriendly=T)
colpal <- "Dark2"
colsfill <- brewer.pal(3,"Dark2")[1:2]
colsAUC <- brewer.pal(10,"RdYlBu")

#Load data for plots----
#All estimates and predictions
dfAll <- data.frame()
#linear
fname <- paste(pathResults,"SimResAppNoteLinear.Rdata",sep="")
load(fname)
dfAll <- rbind(dfAll,df)
#gam
fname <- paste(pathResults,"SimResAppNoteGAM",".Rdata",sep="")
load(fname)
dfAll <- rbind(dfAll,df)
#scam positive
fname <- paste(pathResults,"SimResAppNoteSCAMp",".Rdata",sep="")
load(fname)
dfAll <- rbind(dfAll,df)
#scam positive+monotone increasing
fname <- paste(pathResults,"SimResAppNoteSCAMpmi",".Rdata",sep="")
load(fname)
dfAll <- rbind(dfAll,df)
#hierarchical adaptive discretisation
fname <- paste(pathResults,"SimResAppNoteAD",".Rdata",sep="")
load(fname)
dfAll <- rbind(dfAll,df)

dfAll$G <- as.factor(dfAll$G)
dfAll$method <- factor(dfAll$method, levels=unique(dfAll$method),
                       labels=unique(dfAll$method))
dfAll$Codata <- factor(dfAll$Codata,levels=unique(dfAll$Codata),
                       labels=c("Random","Informative",
                                "Random+transformed","Informative+transformed"))
dfAll$Codata2 <- factor(dfAll$Codata,levels=unique(dfAll$Codata),
                        labels=c("Random","Informative",
                                 "Random","Informative"))

#summarise estimated variance pointwise for continuous co-data per setting
dfEst <- dfAll %>% group_by(G,method,Codata,Covariate,bam.method, Z, ZUntransformed,Transform,Codata2) %>%
  summarise(meanVk = mean(vkfit),q50Vk = quantile(vkfit,0.5),
            q95Vk = quantile(vkfit,0.95),q05Vk = quantile(vkfit,0.05),
            q75Vk = quantile(vkfit,0.75),q25Vk = quantile(vkfit,0.25),
            truevk = mean(truevk)) %>% ungroup()

dfPred1 <- dfAll[dfAll$Covariate==1,] %>% group_by(G,method,Codata,bam.method,Sim,Transform,Codata2) %>% 
  summarise(MSE=mean(MSEridge),method2="ridge")
dfPred2 <- dfAll[dfAll$Covariate==1,] %>% group_by(G,method,Codata,bam.method,Sim,Transform,Codata2) %>% 
  summarise(MSE=mean(MSEecpc),method2="ecpc")
dfPred <- rbind(dfPred1,dfPred2)

#Comparison different smoothing parameter methods
dfGAMs <- data.frame()
#gam
fname <- paste(pathResults,"SimResAppNoteGAM",".Rdata",sep="")
load(fname)
dfGAMs <- rbind(dfGAMs,df)
dfGAMs$G <- as.factor(dfGAMs$G)
dfGAMs$method <- factor(dfGAMs$method, levels=unique(dfGAMs$method),
                        labels=unique(dfGAMs$method))
dfGAMs$Codata <- factor(dfGAMs$Codata,levels=unique(dfGAMs$Codata),
                        labels=c("Random","Informative",
                                 "Random+transformed","Informative+transformed"))
dfGAMs$Codata2 <- factor(dfGAMs$Codata,levels=unique(dfGAMs$Codata),
                         labels=c("Random","Informative",
                                  "Random","Informative"))
#summarise estimated variance pointwise for continuous co-data per setting
dfEstGAMs <- dfGAMs %>% group_by(G,method,Codata,Covariate,bam.method, Z, ZUntransformed,Transform,Codata2) %>%
  summarise(meanVk = mean(vkfit),q50Vk = quantile(vkfit,0.5),
            q95Vk = quantile(vkfit,0.95),q05Vk = quantile(vkfit,0.05),
            q75Vk = quantile(vkfit,0.75),q25Vk = quantile(vkfit,0.25),
            truevk = mean(truevk)) %>% ungroup()

dfPred1GAMs <- dfGAMs[dfGAMs$Covariate==1,] %>% group_by(G,method,Codata,bam.method,Sim,Transform,Codata2) %>% 
  summarise(MSE=mean(MSEridge),method2="ridge")
dfPred2GAMs <- dfGAMs[dfGAMs$Covariate==1,] %>% group_by(G,method,Codata,bam.method,Sim,Transform,Codata2) %>% 
  summarise(MSE=mean(MSEecpc),method2="ecpc")
dfPredGAMs <- rbind(dfPred1GAMs,dfPred2GAMs)


#Table overall mean----
dfTable <- dfPred %>% group_by(G,method,Codata,bam.method,Transform,method2) %>% 
  summarise(MSE=mean(MSE)) %>% ungroup()

#Table average computing times + standard deviation----
timeTable <- dfAll[dfAll$Covariate==1,] %>% group_by(method,bam.method,G) %>% 
  summarise(meanTime=mean(Time),sdTime=sd(Time)) %>% ungroup()
#library("writexl")
#fname <- paste(pathFigures,"TimeTableRaw.xlsx")
#write_xlsx(timeTable,fname)


#Figure estimated variance vs continuous co-data for 1 run, different bam.methods----
g <- G[1]
Sim <- 1

figname <- paste(pathFigures,"Estimates_G",g,"bammethod.pdf",sep="")
pdf(width = wdthpdf*1.2, height = hghtpdf,
    file = figname)
ggplot(dfGAMs[dfGAMs$Sim==Sim & dfGAMs$Transform==F & dfGAMs$G%in%c(1,G) ,])+
  aes(x=Z,col=bam.method)+
  facet_grid(.~Codata,scales="free_x")+
  geom_point(aes(x=Z,y=truevk),alpha=0.1,col="black")+
  geom_line(aes(y=vkfit,linetype=G),size=1.5)+
  labs(y="Prior variance",x="Continuous co-data variable")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
dev.off()

#Figure estimated variance vs continuous co-data for all simulations, fixed G and fixed bam.method----
g <- G[1]
figname <- paste(pathFigures,"Estimates_G",g,"BW2.pdf",sep="")
pdf(width = wdthpdf*1.2, height = hghtpdf*2,
    file = figname)
ggplot(dfEst[dfEst$G%in%c(1,2,g,7) & dfEst$Transform==F &
               dfEst$bam.method%in%c("ML","none","splits"),])+# & dfEst$method=="ML",])+
  aes(x=Z)+
  #facet_wrap(Codata~method,scales="free_x",nrow = 2)+
  facet_grid(method~Codata,scales="free_x")+
  geom_point(aes(x=Z,y=truevk),alpha=0.3,col="grey20")+
  geom_ribbon(aes(ymin=q05Vk,ymax=q95Vk),linetype=0,alpha=0.1)+ 
  geom_ribbon(aes(ymin=q25Vk,ymax=q75Vk),linetype=0,alpha=0.2)+ 
  geom_line(aes(y=q50Vk),size=1)+ #median
  #geom_line(aes(y=meanVk),size=1)+ #mean
  labs(y="Prior variance",x="Continuous co-data variable")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
dev.off()

#Figure estimated variance vs continuous co-data for all simulations, gam only, fixed G and different bam.method----
g <- G[1]
figname <- paste(pathFigures,"Estimates_G",g,"GAM.pdf",sep="")
pdf(width = wdthpdf*1.2, height = hghtpdf*2,
    file = figname)
ggplot(dfEstGAMs[dfEstGAMs$G%in%c(1,2,g) & dfEstGAMs$Transform==F ,])+
  aes(x=Z)+
  facet_grid(bam.method~Codata,scales="free_x")+
  geom_point(aes(x=Z,y=truevk),alpha=0.3,col="black")+
  geom_ribbon(aes(ymin=q05Vk,ymax=q95Vk),linetype=0,alpha=0.1)+ 
  geom_ribbon(aes(ymin=q25Vk,ymax=q75Vk),linetype=0,alpha=0.2)+ 
  geom_line(aes(y=q50Vk),size=1)+ #median
  labs(y="Prior variance",x="Continuous co-data variable")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
dev.off()



#Figure boxplot prediction performance vs co-data setting for all simulations, fixed G and fixed bam.method----
g <- G[1]
temp <- dfPred[dfPred$G%in%c(1,2) & dfPred$bam.method%in%c("none") &
                 dfPred$method2=="ridge" ,]
temp$method <- "ridge"
temp$G <- factor(1,levels=c(1,2,7,20,50),labels=c(1,2,7,20,50))
temp2 <- dfPred[dfPred$G%in%c(1,2,g,7) & dfPred$bam.method%in%c("none","splits","ML") &
                  dfPred$method2=="ecpc" ,]
temp2$G <- factor(temp2$G,levels=c(1,2,7,20,50),labels=c(1,2,7,20,50))
temp2 <- rbind(temp,temp2)
temp2$method <- factor(temp2$method, levels=unique(temp2$method)[c(1,2,4,5,6,3)],
                       labels=unique(temp2$method)[c(1,2,4,5,6,3)])

figname <- paste(pathFigures,"Prediction_G",g,".pdf",sep="")
pdf(width = wdthpdf*1.2, height = hghtpdf,
    file = figname)
ggplot(temp2[temp2$Transform==F,])+#& !(dfPred$Codata=="noninformative"),])+
  aes(x=method,y=MSE)+
  #geom_boxplot(data=temp2[temp2$method=="ridge",],aes(fill=method))+ #ridge prediction performance
  geom_boxplot(fill="grey80")+
  facet_grid(.~Codata2)+
  #coord_cartesian(ylim=c(8,32))+
  #scale_y_log10()+
  labs(y="MSE",x="Co-data model")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts,angle=30,vjust=1,hjust=1),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
dev.off()

#Figure boxplot prediction performance vs co-data setting for all simulations, gam only, fixed G and different bam.method----
g <- G[1]
temp <- dfPred[dfPred$G%in%c(1,2) & dfPred$bam.method%in%c("none") &
                 dfPred$method2=="ridge" ,]
temp$method <- "ridge"
temp$G <- factor(1,levels=c(1,20,50),labels=c(1,20,50))
temp2 <- dfPredGAMs[dfPredGAMs$G%in%c(G) & dfPredGAMs$method%in%c("gam") &
                      dfPredGAMs$method2=="ecpc" ,]
temp2$G <- factor(temp2$G,levels=c(1,20,50),labels=c(1,20,50))
temp2 <- rbind(temp,temp2)
temp2$bam.method <- factor(temp2$bam.method, levels=unique(temp2$bam.method)[c(1,4,2,3,5)],
                           labels=c("ridge",unique(temp2$bam.method)[c(4,2,3,5)]))

figname <- paste(pathFigures,"Prediction_G",g,"GAM.pdf",sep="")
pdf(width = wdthpdf*1.2, height = hghtpdf,
    file = figname)
ggplot(temp2[temp2$Transform==F,])+#& !(dfPred$Codata=="noninformative"),])+
  aes(x=bam.method,y=MSE)+
  #geom_boxplot(data=temp2[temp2$method=="ridge",],aes(fill=method))+ #ridge prediction performance
  geom_boxplot(aes(fill=G))+
  facet_grid(.~Codata2)+
  #scale_y_log10()+-
  labs(y="MSE",x="")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts,angle=30,vjust=1,hjust=1),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
dev.off()




###############################################################################

#Section 3.2 Simulation study on variable selection----
#Load libraries----
library(ecpc)
library(squeezy)
library(glmnet)
library(fwelnet)
library(dplyr) #for data wrangling results
library(ggplot2) #for plotting results
library(RColorBrewer) #for plotting results
library(foreach) #for parallel computing
library(doParallel) #for parallel computing
library(mvtnorm)


#optional: 
if(runParallel){ #set to 1 to setup parallel backend to use many processors
  cores=detectCores()
  if(!("cl"%in%ls())){
    cl <- makeCluster(cores-1) #not to overload your computer
    registerDoParallel(cl)
  }
}

setting_q <- 1
q <- 5/6 #proportion zeros2
maxsel <- c(rep(2:10),10*2:29)
alp_range <- seq(0,1,length.out = length(maxsel))

#Simulate toy data ------
p<-300 #number of covariates
n<-100 #sample size training data set
n2<-100 #sample size test data set
sigma <- 1
tauglobal <- 0.1 #prior variance; 0.1, 0.5
model <- "linear"; fam <- "gaussian"
rho <- 0 #correlation in observed X matrix
Sigma <- matrix(rho,p,p); diag(Sigma) <- 1

#simulate all betas i.i.d. from beta_k~N(mean=0,sd=sqrt(0.1)):
set.seed(7474)
Dat <- list()
Dat$beta <- rnorm(p,sd=sqrt(tauglobal))
ind0 <- sample(1:p,floor(q*p),replace=F)
Dat$beta[ind0] <- 0
Dat$beta <- Dat$beta/sqrt(sum(Dat$beta^2))*sqrt(tauglobal*p)
#Dat$Xctd <- matrix(rnorm(n*p) ,n,p)
#Dat$X2ctd <- matrix(rnorm(n2*p) ,n2,p)
Dat$Xctd <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)
Dat$X2ctd <- rmvnorm(n2,mean=rep(0,p),sigma=Sigma)
Dat$lp2 <- Dat$X2ctd%*%Dat$beta
hist(Dat$lp2)
Dat$Y <- rnorm(n,mean=c(Dat$Xctd%*%Dat$beta),sd=sigma)

#index true, non-zero betas
ind_nonzero <- which(Dat$beta!=0)
ind_zero <- which(Dat$beta==0)

if(generateData){
  AllData <- list()
  for(i in 1:nSim){
    AllData[[i]] <- list()
    AllData[[i]]$beta <- Dat$beta
    if(rho==0){
      AllData[[i]]$Xctd <- matrix(rnorm(n*p) ,n,p)
      AllData[[i]]$X2ctd <- matrix(rnorm(n2*p) ,n2,p)
    }else{
      AllData[[i]]$Xctd <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)
      AllData[[i]]$X2ctd <- rmvnorm(n2,mean=rep(0,p),sigma=Sigma)
    }
    means <- apply(AllData[[i]]$Xctd,2,mean)
    sds <- apply(AllData[[i]]$Xctd,2,sd)
    
    AllData[[i]]$Y <- rnorm(n,mean=c(AllData[[i]]$Xctd%*%AllData[[i]]$beta),sd=sigma)
    AllData[[i]]$Y2 <- rnorm(n2,mean=c(AllData[[i]]$X2ctd%*%AllData[[i]]$beta),sd=sigma)
    
    #standardise data
    #AllData[[i]]$Xctd <- apply(AllData[[i]]$Xctd,2,function(x) (x-mean(x))/sd(x))
    #AllData[[i]]$X2ctd <- t((t(AllData[[i]]$X2ctd)-means)/sds)
  }
  save(Dat,AllData,file=paste("SimDataSparse",setting_q,rho,".Rdata",sep=""))
}
load(paste("SimDataSparse",setting_q,rho,".Rdata",sep=""))

#Co-data settings----
G <- c(20,50) #number of splines
Z.all <- list()
ZI.all <- list() #co-data matrix with intercept
Zs.all <- list()
S1 <- list()
Con.p <- list() #positivity constraints
Con.pmi <- list() #positivity and monotonically increasing constraints
#setting 1: non-informative, unequally spaced
Z.all[["noninformative"]] <- rnorm(p) #for linear co-data model
ZI.all[["noninformative"]] <- cbind(rep(1,p), Z.all[["noninformative"]])
#for generalised additive co-data model:
for(g in 1:length(G)){
  Zs.all[[g]] <- list()
  Zs.all[[g]][["noninformative"]] <- createZforSplines(values=Z.all[["noninformative"]], G=G[g], bdeg=3) 
  S1[[g]] <- createS(orderPen=2, G=G[g]) #create 2nd order difference penalty matrix (same for all co-data)
  
  Con.p[[g]] <- createCon(G=G[g], shape="positive") #create constraints
  Con.pmi[[g]] <- createCon(G=G[g], shape="positive+monotone.i") #create constraints
}
plot(Z.all[[1]],Dat$beta^2)

#setting 2: informative, unequally spaced, information at edge
set.seed(101010)
Z.all[["size.edge"]] <- abs(Dat$beta)+rnorm(p,sd=sd(Dat$beta)/10)
ZI.all[["size.edge"]] <- cbind(rep(1,p), Z.all[["size.edge"]])
#for generalised additive co-data model:
for(g in 1:length(G)){
  Zs.all[[g]][["size.edge"]] <- createZforSplines(values=Z.all[["size.edge"]], G=G[g], bdeg=3) 
}
plot(Z.all[[2]],Dat$beta^2)

#setting 3: informative, unequally spaced, information at 2 edges
set.seed(101010)
Z.all[["size.edge2"]] <- sign(Dat$beta)*abs(Dat$beta)+rnorm(p,sd=sd(Dat$beta)/10)
ZI.all[["size.edge2"]] <- cbind(rep(1,p), Z.all[["size.edge2"]])
#for generalised additive co-data model:
for(g in 1:length(G)){
  Zs.all[[g]][["size.edge2"]] <- createZforSplines(values=Z.all[["size.edge2"]], G=G[g], bdeg=3) 
}
plot(Z.all[[3]],Dat$beta^2)

#fwelnet
Z_fwelnet <- list(matrix(Z.all[["noninformative"]],nrow = p),
                  matrix(Z.all[["size.edge"]],nrow = p),
                  matrix(Z.all[["size.edge2"]],nrow = p))



#Simulation settings----
logname <- "logSimSparse.txt"
write(paste(Sys.time(),logname),file=logname,append=T)


#Fit ecpc+posterior selection on simulated data sets: generalised additive co-data model----
fname <- paste(pathResults,"SimResAppNoteSparse_GAM",setting_q,rho,".Rdata",sep="")
#fname <- paste(pathResults,"SimResAppNoteGAMSparse_splits",q,rho,".Rdata",sep="") #hypershrinkage=ridge
print(fname)

if(run_ecpc){
  #df <- data.frame()
  #for(sim in 1:nSim){
  finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
                         .packages = c("ecpc", "squeezy")) %dopar% {
                           df <- data.frame()
                           dfSelect <- data.frame()
                           for(setting in 1:3){
                             for(g in 1){
                               for(method in c("ML")){
                                 tic<-proc.time()[[3]]
                                 
                                 if(method=="splits"){
                                   fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
                                               Z=Zs.all[[g]][setting],
                                               paraPen=list(Z1=list(S1=S1[[g]])),
                                               intrcpt.bam = F,
                                               model=model,maxsel=maxsel,
                                               Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
                                               hypershrinkage = "ridge",
                                               est_beta_method = "multiridge")
                                 }else{
                                   fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
                                               Z=Zs.all[[g]][setting],
                                               paraPen=list(Z1=list(S1=S1[[g]])),
                                               bam.method=method,intrcpt.bam = F,
                                               model=model,maxsel=maxsel,
                                               Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
                                               est_beta_method = "multiridge")
                                 }
                                 
                                 
                                 toc <- proc.time()[[3]]-tic
                                 
                                 vk <- (Zs.all[[g]][[setting]]%*%fit$gamma )*fit$tauglobal
                                 #plot(Z.all[[setting]],vk)
                                 
                                 temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
                                 temp$Z <- Z.all[[setting]]
                                 temp$Time <- toc
                                 temp$Covariate <- 1:p
                                 temp$G <- G[g]
                                 temp$setting <- setting
                                 temp$bam.method <- method
                                 temp$method <- "gam"
                                 temp$MSEridge <- fit$MSEridge
                                 temp$MSEecpc <- fit$MSEecpc
                                 temp$Sim <- sim
                                 temp$Codata <- names(Z.all)[setting]
                                 temp$Transform <- F; if(setting>=3) temp$Transform <- T
                                 temp$ZUntransformed <- Z.all[[setting]]
                                 temp$q <- q
                                 temp$rho <- rho
                                 if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
                                 df <- rbind(df,temp)
                                 
                                 #sensitivity = true positive rate = true positives/total positives
                                 sensitivity <- apply(fit$betaPost,2,function(x){
                                   ind_nonzero_est <- which(x!=0)
                                   TPR <-sum(ind_nonzero_est%in%ind_nonzero)/length(ind_nonzero)
                                   return(TPR)
                                 })
                                 #precision = true negative rate = true negatives/total negatives
                                 precision <- apply(fit$betaPost,2,function(x){
                                   ind_zero_est <- which(x==0)
                                   TNR <-sum(ind_zero_est%in%ind_zero)/length(ind_zero)
                                   return(TNR)
                                 })
                                 lp_train <- AllData[[sim]]$Xctd%*%fit$betaPost
                                 MSEtrain <- apply(lp_train,2,function(lp) mean((lp-AllData[[sim]]$Y)^2))
                                 
                                 temp<-data.frame("Tuningparam"=maxsel, sensitivity, precision)
                                 temp$TypeTuning <- "#params"
                                 temp$Time <- toc
                                 temp$G <- G[g]
                                 temp$setting <- setting
                                 temp$bam.method <- method
                                 temp$method <- "gam"
                                 temp$MSEridge <- fit$MSEridge
                                 temp$MSEecpc <- c(fit$MSEPost)
                                 temp$MSEtrain <- c(MSEtrain)
                                 if(is.null(fit$MSEPost)) temp$MSEecpc <- rep(NaN,length(maxsel))
                                 temp$Sim <- sim
                                 temp$Codata <- names(Z.all)[setting]
                                 temp$Transform <- F; if(setting>=3) temp$Transform <- T
                                 temp$q <- q
                                 temp$setting_q <- setting_q
                                 temp$rho <- rho
                                 dfSelect <- rbind(dfSelect, temp)
                                 
                                 #fit squeezy
                                 #Use squeezy function to transform estimated ridge penalties to elastic net penalties
                                 notinf <- fit$penalties!=Inf
                                 fit.EN <- lapply(alp_range,function(alp){
                                   fit.EN <- squeezy(Y=AllData[[sim]]$Y, X=AllData[[sim]]$Xctd[,notinf], 
                                                     groupset=lapply(1:sum(notinf), function(x) x), 
                                                     alpha=alp, 
                                                     Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd[,notinf],
                                                     lambdas=fit$penalties[notinf],
                                                     sigmasq=fit$sigmahat)})
                                 beta_squeezy <- sapply(1:length(alp_range),function(i){
                                   betas <- rep(0,p)
                                   betas[notinf] <- fit.EN[[i]]$betaApprox
                                   return(betas)
                                 })
                                 MSE_squeezy <- sapply(1:length(alp_range),function(i){
                                   fit.EN[[i]]$MSEApprox
                                 })
                                 lp_train <- AllData[[sim]]$Xctd%*%beta_squeezy
                                 MSEtrain <- apply(lp_train,2,function(lp) mean((lp-AllData[[sim]]$Y)^2))
                                 
                                 #sensitivity = true positive rate = true positives/total positives
                                 sensitivity <- apply(beta_squeezy,2,function(x){
                                   ind_nonzero_est <- which(x!=0)
                                   TPR <-sum(ind_nonzero_est%in%ind_nonzero)/length(ind_nonzero)
                                   return(TPR)
                                 })
                                 #precision = true negative rate = true negatives/total negatives
                                 precision <- apply(beta_squeezy,2,function(x){
                                   ind_zero_est <- which(x==0)
                                   TNR <-sum(ind_zero_est%in%ind_zero)/length(ind_zero)
                                   return(TNR)
                                 })
                                 
                                 temp<-data.frame("Tuningparam"=alp_range, sensitivity, precision)
                                 temp$TypeTuning <- "alpha"
                                 temp$Time <- toc
                                 temp$G <- G[g]
                                 temp$setting <- setting
                                 temp$bam.method <- method
                                 temp$method <- "gam"
                                 temp$MSEridge <- fit$MSEridge
                                 temp$MSEecpc <- MSE_squeezy
                                 temp$MSEtrain <- c(MSEtrain)
                                 temp$Sim <- sim
                                 temp$Codata <- names(Z.all)[setting]
                                 temp$Transform <- F; if(setting>=3) temp$Transform <- T
                                 temp$q <- q
                                 temp$setting_q <- setting_q
                                 temp$rho <- rho
                                 dfSelect <- rbind(dfSelect,temp)
                                 
                               }
                             }
                           }
                           write(paste(Sys.time(),"simulation",sim,"of",nSim,"done"),file=logname,append=T)
                           list("df"=df,"dfSelect"=dfSelect)
                         }                                     
  #str(finalMatrix)
  df2 <- lapply(1:nSim,function(i) finalMatrix[i,1][[1]])
  dfSelect2 <- lapply(1:nSim,function(i) finalMatrix[i,2][[1]])
  df <- df2[[1]]; for(i in 2:nSim) df <- rbind(df,df2[[i]])
  dfSelect <- dfSelect2[[1]]; for(i in 2:nSim) dfSelect <- rbind(dfSelect,dfSelect2[[i]])
  save(finalMatrix,df,dfSelect,file=fname)
}


#Fit glmnet on simulated data sets----
fname <- paste(pathResults,"SimResAppNoteSparse_glmnet",setting_q,rho,".Rdata",sep="")
#fname <- paste(pathResults,"SimResAppNoteGAMSparse_splits",q,rho,".Rdata",sep="") #hypershrinkage=ridge
print(fname)

if(run_glmnet){
  #df <- data.frame()
  #for(sim in 1:nSim){
  finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
                         .packages = c("glmnet")) %dopar% {
                           dfSelect <- data.frame()
                           for(setting in 1:3){
                             for(alp in alp_range){
                               tic<-proc.time()[[3]]
                               
                               fit.glmnet <- glmnet::cv.glmnet(y=AllData[[sim]]$Y,x=AllData[[sim]]$Xctd,
                                                               family=fam,alpha=alp)
                               beta.glmnet <- coef(fit.glmnet,s="lambda.min", exact=TRUE)
                               Ypred.glmnet <- c(predict(fit.glmnet, newx = AllData[[sim]]$X2ctd, 
                                                         s = "lambda.min", type="response", exact=TRUE))
                               MSE.glmnet <- mean((Ypred.glmnet-AllData[[sim]]$Y2)^2)
                               lp_train <- AllData[[sim]]$Xctd%*%beta.glmnet[-1] + beta.glmnet[1] 
                               MSEtrain <- mean((lp_train-AllData[[sim]]$Y)^2)
                               toc <- proc.time()[[3]]-tic
                               
                               #sensitivity = true positive rate = true positives/total positives
                               ind_nonzero_est <- which(beta.glmnet[-1]!=0)
                               sensitivity <- sum(ind_nonzero_est%in%ind_nonzero)/length(ind_nonzero)
                               #precision = true negative rate = true negatives/total negatives
                               ind_zero_est <- which(beta.glmnet[-1]==0)
                               precision <- sum(ind_zero_est%in%ind_zero)/length(ind_zero)
                               
                               temp<-data.frame("Tuningparam"=alp, sensitivity, precision)
                               temp$TypeTuning <- "alpha"
                               temp$Time <- toc
                               temp$G <- 1
                               temp$setting <- setting
                               temp$bam.method <- "none"
                               temp$method <- "glmnet"
                               temp$MSEridge <- NaN
                               temp$MSEecpc <- MSE.glmnet
                               temp$MSEtrain <- MSEtrain
                               temp$Sim <- sim
                               temp$Codata <- names(Z.all)[setting]
                               temp$Transform <- F; if(setting>=3) temp$Transform <- T
                               temp$q <- q
                               temp$setting_q <- setting_q
                               temp$rho <- rho
                               dfSelect <- rbind(dfSelect,temp)
                               
                             }
                           }
                           write(paste(Sys.time(),"simulation",sim,"of",nSim,"done"),file=logname,append=T)
                           list("dfSelect"=dfSelect)
                         }                                     
  #str(finalMatrix)
  dfSelect <- finalMatrix[[1]]; for(i in 2:nSim) dfSelect <- rbind(dfSelect,finalMatrix[[i]])
  save(finalMatrix,dfSelect,file=fname)
}

#Fit fwelnet on simulated data sets----
fname <- paste(pathResults,"SimResAppNoteSparse_fwelnet",setting_q,rho,".Rdata",sep="")
#fname <- paste(pathResults,"SimResAppNoteGAMSparse_splits",q,rho,".Rdata",sep="") #hypershrinkage=ridge
print(fname)

if(run_fwelnet){
  #df <- data.frame()
  #for(sim in 1:nSim){
  finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
                         .packages = c("fwelnet")) %dopar% {
                           dfSelect <- data.frame()
                           for(setting in 1:3){
                             for(alp in alp_range){
                               tic<-proc.time()[[3]]
                               
                               res.fwEN <-  cv.fwelnet(y=AllData[[sim]]$Y,x=AllData[[sim]]$Xctd,
                                                       z=Z_fwelnet[[setting]],
                                                       alpha=alp,family=fam,standardize=F)
                               
                               ind.minlam <- which(res.fwEN$lambda==res.fwEN$lambda.min)
                               
                               betafwEN <- res.fwEN$glmfit$beta[,ind.minlam]
                               a0fwEN <- res.fwEN$glmfit$a0[ind.minlam]
                               Ypred_fwelnet <- AllData[[sim]]$X2ctd%*%betafwEN+a0fwEN
                               MSE_fwelnet <- mean((Ypred_fwelnet-AllData[[sim]]$Y2)^2)
                               lp_train <- AllData[[sim]]$Xctd%*%betafwEN + a0fwEN 
                               MSEtrain <- mean((lp_train-AllData[[sim]]$Y)^2)
                               
                               toc <- proc.time()[[3]]-tic
                               
                               #sensitivity = true positive rate = true positives/total positives
                               ind_nonzero_est <- which(betafwEN!=0)
                               sensitivity <- sum(ind_nonzero_est%in%ind_nonzero)/length(ind_nonzero)
                               #precision = true negative rate = true negatives/total negatives
                               ind_zero_est <- which(betafwEN==0)
                               precision <- sum(ind_zero_est%in%ind_zero)/length(ind_zero)
                               
                               temp<-data.frame("Tuningparam"=alp, sensitivity, precision)
                               temp$TypeTuning <- "alpha"
                               temp$Time <- toc
                               temp$G <- 1
                               temp$setting <- setting
                               temp$bam.method <- "none"
                               temp$method <- "fwelnet"
                               temp$MSEridge <- NaN
                               temp$MSEecpc <- MSE_fwelnet
                               temp$MSEtrain <- c(MSEtrain)
                               temp$Sim <- sim
                               temp$Codata <- names(Z.all)[setting]
                               temp$Transform <- F; if(setting>=3) temp$Transform <- T
                               temp$q <- q
                               temp$setting_q <- setting_q
                               temp$rho <- rho
                               dfSelect <- rbind(dfSelect,temp)
                               
                             }
                           }
                           write(paste(Sys.time(),"simulation",sim,"of",nSim,"done"),file=logname,append=T)
                           list("dfSelect"=dfSelect)
                         }                                     
  #str(finalMatrix)
  dfSelect <- finalMatrix[[1]]; for(i in 2:nSim) dfSelect <- rbind(dfSelect,finalMatrix[[i]])
  save(finalMatrix,dfSelect,file=fname)
}

#Plots: general parameters----
wdth<-600
hght<-wdth*5/8
wdthpdf <- wdth/75
hghtpdf <- hght/75
ts <- 16 #basis text size in figures
ls <- 1.5 #basis line size in figures
ps <- 2 #basis point size in figures
sz <- 2 #point size
strk <- 1.5 #stroke size
palette <- "Dark2"
#display.brewer.all(10,colorblindFriendly=T)
colpal <- "Dark2"
colsfill <- brewer.pal(3,"Dark2")[1:2]
colsAUC <- brewer.pal(10,"RdYlBu")

#Load data for plots----
#All estimates and predictions
dfAll <- data.frame()
#gam
fname <- paste(pathResults,"SimResAppNoteSparse_GAM",setting_q,rho,".Rdata",sep="")
load(fname)
dfSelect$method[dfSelect$TypeTuning=="#params"] <- "ecpc+postselection"
dfSelect$method[dfSelect$TypeTuning=="alpha"] <- "ecpc+squeezy"
dfAll <- rbind(dfAll,dfSelect)
#glmnet
fname <- paste(pathResults,"SimResAppNoteSparse_glmnet",setting_q,rho,".Rdata",sep="")
load(fname)
dfAll <- rbind(dfAll,dfSelect)

#fwelnet
fname <- paste(pathResults,"SimResAppNoteSparse_fwelnet",setting_q,rho,".Rdata",sep="")
load(fname)
dfAll <- rbind(dfAll,dfSelect)

dfAll$Tuningparam2 <- dfAll$Tuningparam
dfAll$TypeTuning2 <- dfAll$TypeTuning
dfAll$Tuningparam2[dfAll$TypeTuning=="#params"] <- (p-dfAll$Tuningparam2[dfAll$TypeTuning=="#params"])/p
dfAll$TypeTuning2[dfAll$TypeTuning=="#params"] <- "proportion zeros" 

dfAll$G <- as.factor(dfAll$G)
dfAll$method <- factor(dfAll$method, levels=unique(dfAll$method)[c(3,4,2,1)],
                       labels=unique(dfAll$method)[c(3,4,2,1)])
dfAll$Codata <- factor(dfAll$Codata,levels=unique(dfAll$Codata),
                       labels=c("Random","Informative+monotone","Informative+convex"))

dfPred <- dfAll %>% group_by(method,Codata,TypeTuning,Tuningparam,TypeTuning2,Tuningparam2) %>% 
  summarise(meanMSE=mean(MSEecpc),method2="EN",
            q50MSE = quantile(MSEecpc,0.5),
            q95MSE = quantile(MSEecpc,0.95),q05MSE = quantile(MSEecpc,0.05),
            q75MSE = quantile(MSEecpc,0.75),q25MSE = quantile(MSEecpc,0.25),
            meanMSEtrain=mean(MSEtrain),method2="EN",
            q50MSEtrain = quantile(MSEtrain,0.5),
            q95MSEtrain = quantile(MSEtrain,0.95),q05MSEtrain = quantile(MSEtrain,0.05),
            q75MSEtrain = quantile(MSEtrain,0.75),q25MSEtrain = quantile(MSEtrain,0.25),
            meansensitivity=mean(sensitivity),
            q50sens = quantile(sensitivity,0.5,na.rm=T),
            q95sens = quantile(sensitivity,0.95,na.rm=T),q05sens = quantile(sensitivity,0.05,na.rm=T),
            q75sens = quantile(sensitivity,0.75,na.rm=T),q25sens = quantile(sensitivity,0.25,na.rm=T),
            meanprecision=mean(precision),
            q50prec = quantile(precision,0.5,na.rm=T),
            q95prec = quantile(precision,0.95,na.rm=T),q05prec = quantile(precision,0.05,na.rm=T),
            q75prec = quantile(precision,0.75,na.rm=T),q25prec = quantile(precision,0.25,na.rm=T)) %>% ungroup()


#Figure MSE for all simulations, tuning parameter 2----
figname <- paste(pathFigures,"Prediction_sparseMSE.pdf",sep="")
pdf(width = wdthpdf*1.8, height = hghtpdf,
    file = figname)
ggplot(dfPred)+# & dfEst$method=="ML",])+
  aes(x=Tuningparam2, col=method, fill=method)+
  facet_grid(.~Codata)+
  geom_ribbon(aes(ymin=q05MSE,ymax=q95MSE),linetype=0,alpha=0.1)+ 
  geom_ribbon(aes(ymin=q25MSE,ymax=q75MSE),linetype=0,alpha=0.2)+ 
  #geom_line(aes(y=q50MSE, linetype=TypeTuning2),size=1)+ #median
  geom_line(aes(y=meanMSE, linetype=TypeTuning2),size=1)+ #mean
  #scale_y_log10()+
  labs(y="MSE",x="Tuning parameter")+
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"),
         linetype=guide_legend(title="Type tuning parameter"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
dev.off()

#Figure sensitivity for all simulations, tuning parameter 2----
figname <- paste(pathFigures,"Prediction_sparseSens.pdf",sep="")
pdf(width = wdthpdf*1.2, height = hghtpdf*2,
    file = figname)
ggplot(dfPred)+# & dfEst$method=="ML",])+
  aes(x=Tuningparam2, linetype=TypeTuning2, col=method, fill=method)+
  facet_grid(.~Codata)+
  #geom_line(aes(y=MSE),alpha=1, size=ls)+ #mean
  geom_ribbon(aes(ymin=q05sens,ymax=q95sens),linetype=0,alpha=0.1)+ 
  geom_ribbon(aes(ymin=q25sens,ymax=q75sens),linetype=0,alpha=0.2)+ 
  geom_line(aes(y=q50sens),size=1)+ #median
  #geom_line(aes(y=meansensitivity),size=1)+ #mean
  labs(y="Sensitivity",x="Tuning parameter")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
dev.off()


#Figure precision for all simulations----
figname <- paste(pathFigures,"Prediction_sparsePrec.pdf",sep="")
pdf(width = wdthpdf*1.2, height = hghtpdf*2,
    file = figname)
ggplot(dfPred)+# & dfEst$method=="ML",])+
  aes(x=Tuningparam2, linetype=TypeTuning2, col=method, fill=method)+
  facet_grid(.~Codata)+
  #geom_line(aes(y=MSE),alpha=1, size=ls)+ #mean
  geom_ribbon(aes(ymin=q05prec,ymax=q95prec),linetype=0,alpha=0.1)+ 
  geom_ribbon(aes(ymin=q25prec,ymax=q75prec),linetype=0,alpha=0.2)+ 
  geom_line(aes(y=q50prec),size=1)+ #median
  #geom_line(aes(y=meanVk),size=1)+ #mean
  labs(y="Precision",x="Tuning parameter")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
dev.off()

#Figure sensitivity vs precision for all simulations----
figname <- paste(pathFigures,"Prediction_sparseSensPrec.pdf",sep="")
pdf(width = wdthpdf*1.8, height = hghtpdf,
    file = figname)
ggplot(dfPred)+# & dfEst$method=="ML",])+
  aes(col=method, shape=TypeTuning2)+
  facet_grid(.~Codata)+
  geom_line(aes(y=meansensitivity,x=1-meanprecision),size=ls, alpha=0.2)+ #mean
  geom_point(aes(y=meansensitivity,x=1-meanprecision),size=ps,alpha=0.8)+ #mean
  guides(color=guide_legend(title="Method"),
         shape=guide_legend(title="Type tuning parameter"))+
  #geom_line(aes(y=meanVk),size=1)+ #mean
  labs(y="Sensitivity",x="1-precision")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
dev.off()


###############################################################################

#Section 4: Application to lymph node metastasis classification################
################################################################################

#Load packages----
library(ecpc)
library(devtools)
library(CoRF) #load package
library(ggplot2)
library(dplyr)
library(pROC)
library(ggpubr)

#Load data from CoRF package----
data(LNM_Example)
#make response of RespValidation numeric such that its format matches that of RespTrain
RespValidationNum <- as.numeric(RespValidation)-1 
maxSel <- c(2:10,10*(2:10)) #maximum selected variables
G_all <- c(20,50) #use 20 or 50 spline basis functions
n <- dim(TrainData)[1]
p <- dim(TrainData)[2]

#Example in paper----


if(run_example){
  set.seed(3)
  G <- G_all[1] #number of splines
  #co-data source 1: genes corresponding to a known signature or not
  GroupsetSig <- createGroupset(as.factor(CoDataTrain$RoepmanGenes)) #list of groups
  Z_sig <- createZforGroupset(GroupsetSig) #co-data matrix
  colnames(Z_sig) <- c("out","in")
  
  #Co-data 2: correlation
  Z_cor <- createZforSplines(values=CoDataTrain$Corrs, G=G)
  S1_cor <- createS(orderPen=2, G=G) #create difference penalty matrix
  Con_cor <- createCon(G=G, shape="positive+monotone.i")
  
  #Co-data 3: p-values of external study using a different measurement technique
  #2 probes have NA p-value, set to maximum observed p-value (approx. 1)
  CoDataTrain$pvalsVUmc[is.na(CoDataTrain$pvalsVUmc)] <- max(CoDataTrain$pvalsVUmc,na.rm=TRUE)
  Z_pvals <- createZforSplines(values=CoDataTrain$pvalsVUmc, G=G)
  S1_pvals <- createS(G=G) #create difference penalty matrix
  Con_pvals <- createCon(G=G, shape="positive+monotone.d")
  
  
  #save continuous co-data values in seperate list for plotting
  values <- list("Signature" = NULL, 
                 "Correlation" = CoDataTrain$Corrs,
                 "P-values" = CoDataTrain$pvalsVUmc)
  
  Z_all <- list("Signature" = Z_sig, 
                "Correlation" = Z_cor,
                "P-values" = Z_pvals)
  
  #Fit model: all three co-data matrices, with constraints
  fname <- "Res.Rdata"
  
  tic <- proc.time()[[3]]
  Res<-ecpc(Y=RespTrain, X=TrainData, #training data
            Z=Z_all, 
            paraPen = list(Z2=list(S1=S1_cor), Z3=list(S1=S1_pvals)),
            paraCon = list(Z2=Con_cor, Z3=Con_pvals),
            Y2=RespValidationNum, X2=ValidationData, #test data
            maxsel=maxSel) #maximum number of selected covariates for posterior selection
  toc <- proc.time()[[3]]; 
  print(toc - tic) #elapsed time
  Res$time <- toc-tic
  
  #visualise estimates
  plot(Res, show = "priorweights", Z=Z_all, values=values)
  
  #obtain predictions (same as Res$Ypred)
  Ypred <- predict(Res, X2=ValidationData)
  
  #Fit posterior selected models (same as Res$betaPost)
  sparseModels <- postSelect(Res, X=TrainData, Y=RespTrain, maxsel=maxSel)
  
  
  #concatenate results on predictions in data frame
  ntest <- length(RespValidation) #number of test samples in the validation data
  p<-dim(CoDataTrain)[1] #number of genes
  df<-data.frame("Ypred"= c(Res$Ypred, Res$Ypredridge, c(Res$YpredPost)),
                 "Truth" = rep(RespValidationNum, 2+length(maxSel)),
                 "Sample" = rep(1:ntest,2+length(maxSel)),
                 "NumberSelectedVars" = rep(c(sum(Res$beta!=0),p,maxSel),each=ntest),
                 "Method" = as.factor(rep(c("ecpc","ordinary.ridge",rep("ecpc",length(maxSel))),each=ntest)),
                 "Sparsity" = as.factor(rep(c("dense","dense", rep("sparse",length(maxSel))),each=ntest)),
                 "MethodxNumberVars" = as.factor(rep(c("ecpc","ordinary.ridge",
                                                       paste("ecpc",maxSel,"vars",sep="")),each=ntest)))
  
  #Compute AUC for each method and store in a data frame
  dfAUC <- data.frame()
  for (i in levels(df$MethodxNumberVars)) {
    temp <- data.frame("Method" = unique(df$Method[df$MethodxNumberVars==i]),
                       "Sparsity"=unique(df$Sparsity[df$MethodxNumberVars==i]),
                       "MethodxNumberVars"=i)
    rocpROC <- pROC::roc(predictor = df$Ypred[df$MethodxNumberVars == i], 
                         response = df$Truth[df$MethodxNumberVars == i], 
                         smooth = F, auc = T, levels = c(0, 1), direction = "<")
    temp$AUC <- rocpROC$auc[1]
    temp$NumberSelectedVars <- mean(df$NumberSelectedVars[df$MethodxNumberVars == i])
    dfAUC <- rbind(dfAUC, temp)
  }
  lims <- range(dfAUC$AUC,0.5) #set limits for AUC plots
  dfAUC$Codata <- "all"
  
  save(Res,df,dfAUC,Z_all,values,file=fname)
}else{
  #example
  fname <- "Res.Rdata"
  load(fname)
}
Res_example <- Res
Z_example <- Z_all
values_example <- values


#Run all analyses and save results for plotting----

Z_all <- list(); values <- list()
Z_all_rank <- list(); values_rank <- list()
for(g in 1:length(G_all)){
  #Format co-data----
  G <- G_all[g] #number of splines
  #co-data source 1: genes corresponding to a known signature or not
  GroupsetSig <- createGroupset(as.factor(CoDataTrain$RoepmanGenes)) #list of groups
  Z_sig <- createZforGroupset(GroupsetSig) #co-data matrix
  colnames(Z_sig) <- c("out","in")
  
  #Co-data 2: correlation
  Z_cor <- createZforSplines(values=CoDataTrain$Corrs, G=G)
  S1_cor <- createS(orderPen=2, G=G) #create difference penalty matrix
  Con_cor <- createCon(G=G, shape="positive+monotone.i")
  Con_cor2 <- createCon(G=G, shape="positive")
  
  #Co-data 2a: rank correlation
  ecdfZ <- ecdf(CoDataTrain$Corrs)
  ZRankCorrs <- ecdfZ(CoDataTrain$Corrs)
  Z_cor_rank <- createZforSplines(values=ZRankCorrs, G=G)
  S1_cor_rank <- createS(G=G) #create difference penalty matrix
  
  #Co-data 3: p-values of external study using a different measurement technique
  #2 probes have NA p-value, set to maximum observed p-value (approx. 1)
  CoDataTrain$pvalsVUmc[is.na(CoDataTrain$pvalsVUmc)] <- max(CoDataTrain$pvalsVUmc,na.rm=TRUE)
  Z_pvals <- createZforSplines(values=CoDataTrain$pvalsVUmc, G=G)
  S1_pvals <- createS(G=G) #create difference penalty matrix
  Con_pvals <- createCon(G=G, shape="positive+monotone.d")
  Con_pvals2 <- createCon(G=G, shape="positive")
  
  #Co-data 3a: rank p-values
  ZRankPvals <- ecdf(CoDataTrain$pvalsVUmc)(CoDataTrain$pvalsVUmc)
  Z_pvals_rank <- createZforSplines(values=ZRankPvals, G=G)
  S1_pvals_rank <- createS(G=G) #create difference penalty matrix
  
  
  #save continuous co-data values in seperate list for plotting
  values[[g]] <- list("Signature" = NULL, 
                      "Correlation" = CoDataTrain$Corrs,
                      "P-values" = CoDataTrain$pvalsVUmc)
  values_rank[[g]] <- list("Signature" = NULL, 
                           "Correlation" = ZRankCorrs,
                           "P-values" = ZRankPvals)
  
  Z_all[[g]] <- list("Signature" = Z_sig, 
                     "Correlation" = Z_cor,
                     "P-values" = Z_pvals)
  
  Z_all_rank[[g]] <- list("Signature" = Z_sig, 
                          "Correlation" = Z_cor_rank,
                          "P-values" = Z_pvals_rank)
  
  
  
  #Fit model 1: all three co-data matrices, with constraints----
  fname <- paste("Res1G",g,".Rdata",sep="")
  print(fname)
  
  if(run_withConstraints){
    tic <- proc.time()[[3]]
    Res<-ecpc(Y=RespTrain, X=TrainData, #training data
              Z=Z_all[[g]], 
              paraPen = list(Z2=list(S1=S1_cor), Z3=list(S1=S1_pvals)),
              paraCon = list(Z2=Con_cor, Z3=Con_pvals),
              Y2=RespValidationNum, X2=ValidationData, #test data
              maxsel=maxSel) #maximum number of selected covariates for posterior selection
    toc <- proc.time()[[3]]; 
    print(toc - tic) #elapsed time
    Res$time <- toc-tic
    
    plot(Res, show = "priorweights", Z=Z_all[[g]], values=values[[g]])
    
    #concatenate results on predictions in data frame
    ntest <- length(RespValidation) #number of test samples in the validation data
    p<-dim(CoDataTrain)[1] #number of genes
    df<-data.frame("Ypred"= c(Res$Ypred, Res$Ypredridge, c(Res$YpredPost)),
                   "Truth" = rep(RespValidationNum, 2+length(maxSel)),
                   "Sample" = rep(1:ntest,2+length(maxSel)),
                   "NumberSelectedVars" = rep(c(sum(Res$beta!=0),p,maxSel),each=ntest),
                   "Method" = as.factor(rep(c("ecpc","ordinary.ridge",rep("ecpc",length(maxSel))),each=ntest)),
                   "Sparsity" = as.factor(rep(c("dense","dense", rep("sparse",length(maxSel))),each=ntest)),
                   "MethodxNumberVars" = as.factor(rep(c("ecpc","ordinary.ridge",
                                                         paste("ecpc",maxSel,"vars",sep="")),each=ntest)))
    
    #Compute AUC for each method and store in a data frame
    dfAUC <- data.frame()
    for (i in levels(df$MethodxNumberVars)) {
      temp <- data.frame("Method" = unique(df$Method[df$MethodxNumberVars==i]),
                         "Sparsity"=unique(df$Sparsity[df$MethodxNumberVars==i]),
                         "MethodxNumberVars"=i)
      rocpROC <- pROC::roc(predictor = df$Ypred[df$MethodxNumberVars == i], 
                           response = df$Truth[df$MethodxNumberVars == i], 
                           smooth = F, auc = T, levels = c(0, 1), direction = "<")
      temp$AUC <- rocpROC$auc[1]
      temp$NumberSelectedVars <- mean(df$NumberSelectedVars[df$MethodxNumberVars == i])
      dfAUC <- rbind(dfAUC, temp)
    }
    lims <- range(dfAUC$AUC,0.5) #set limits for AUC plots
    dfAUC$Codata <- "all"
    
    save(Res,df,dfAUC,file=fname)
  }
  
  
  #Fit model 2: all three co-data matrices, no constraints----
  fname <- paste("Res2G",g,".Rdata",sep="")
  print(fname)
  
  if(run_withoutConstraints){
    tic <- proc.time()[[3]]
    Res<-ecpc(Y=RespTrain, X=TrainData, #training data
              Z=Z_all[[g]], 
              paraPen = list(Z2=list(S1=S1_cor), Z3=list(S1=S1_pvals)),
              Y2=RespValidationNum, X2=ValidationData, #test data
              maxsel=maxSel) #maximum number of selected covariates for posterior selection
    toc <- proc.time()[[3]]; 
    print(toc - tic) #elapsed time
    Res$time <- toc-tic
    
    plot(Res, show = "priorweights", Z=Z_all[[g]], values=values[[g]])
    
    #concatenate results on predictions in data frame
    ntest <- length(RespValidation) #number of test samples in the validation data
    p<-dim(CoDataTrain)[1] #number of genes
    df<-data.frame("Ypred"= c(Res$Ypred, Res$Ypredridge, c(Res$YpredPost)),
                   "Truth" = rep(RespValidationNum, 2+length(maxSel)),
                   "Sample" = rep(1:ntest,2+length(maxSel)),
                   "NumberSelectedVars" = rep(c(sum(Res$beta!=0),p,maxSel),each=ntest),
                   "Method" = as.factor(rep(c("ecpc","ordinary.ridge",rep("ecpc",length(maxSel))),each=ntest)),
                   "Sparsity" = as.factor(rep(c("dense","dense", rep("sparse",length(maxSel))),each=ntest)),
                   "MethodxNumberVars" = as.factor(rep(c("ecpc","ordinary.ridge",
                                                         paste("ecpc",maxSel,"vars",sep="")),each=ntest)))
    
    #Compute AUC for each method and store in a data frame
    dfAUC <- data.frame()
    for (i in levels(df$MethodxNumberVars)) {
      temp <- data.frame("Method" = unique(df$Method[df$MethodxNumberVars==i]),
                         "Sparsity"=unique(df$Sparsity[df$MethodxNumberVars==i]),
                         "MethodxNumberVars"=i)
      rocpROC <- pROC::roc(predictor = df$Ypred[df$MethodxNumberVars == i], 
                           response = df$Truth[df$MethodxNumberVars == i], 
                           smooth = F, auc = T, levels = c(0, 1), direction = "<")
      temp$AUC <- rocpROC$auc[1]
      temp$NumberSelectedVars <- mean(df$NumberSelectedVars[df$MethodxNumberVars == i])
      dfAUC <- rbind(dfAUC, temp)
    }
    lims <- range(dfAUC$AUC,0.5) #set limits for AUC plots
    dfAUC$Codata <- "all"
    
    save(Res,df,dfAUC,file=fname)
  }
  
  
  #Fit model 3: all three co-data matrices, no constraints, with transformed scale----
  fname <- paste("Res3G",g,".Rdata",sep="")
  print(fname)
  
  if(run_transformed){
    tic <- proc.time()[[3]]
    Res<-ecpc(Y=RespTrain, X=TrainData, #training data
              Z=Z_all_rank[[g]], 
              paraPen = list(Z2=list(S1=S1_cor), Z3=list(S1=S1_pvals)),
              Y2=RespValidationNum, X2=ValidationData, #test data
              maxsel=maxSel) #maximum number of selected covariates for posterior selection
    toc <- proc.time()[[3]]; 
    print(toc - tic) #elapsed time
    Res$time <- toc-tic
    
    plot(Res, show = "priorweights", Z=Z_all_rank[[g]], values=values[[g]])
    
    #concatenate results on predictions in data frame
    ntest <- length(RespValidation) #number of test samples in the validation data
    p<-dim(CoDataTrain)[1] #number of genes
    df<-data.frame("Ypred"= c(Res$Ypred, Res$Ypredridge, c(Res$YpredPost)),
                   "Truth" = rep(RespValidationNum, 2+length(maxSel)),
                   "Sample" = rep(1:ntest,2+length(maxSel)),
                   "NumberSelectedVars" = rep(c(sum(Res$beta!=0),p,maxSel),each=ntest),
                   "Method" = as.factor(rep(c("ecpc","ordinary.ridge",rep("ecpc",length(maxSel))),each=ntest)),
                   "Sparsity" = as.factor(rep(c("dense","dense", rep("sparse",length(maxSel))),each=ntest)),
                   "MethodxNumberVars" = as.factor(rep(c("ecpc","ordinary.ridge",
                                                         paste("ecpc",maxSel,"vars",sep="")),each=ntest)))
    
    #Compute AUC for each method and store in a data frame
    dfAUC <- data.frame()
    for (i in levels(df$MethodxNumberVars)) {
      temp <- data.frame("Method" = unique(df$Method[df$MethodxNumberVars==i]),
                         "Sparsity"=unique(df$Sparsity[df$MethodxNumberVars==i]),
                         "MethodxNumberVars"=i)
      rocpROC <- pROC::roc(predictor = df$Ypred[df$MethodxNumberVars == i], 
                           response = df$Truth[df$MethodxNumberVars == i], 
                           smooth = F, auc = T, levels = c(0, 1), direction = "<")
      temp$AUC <- rocpROC$auc[1]
      temp$NumberSelectedVars <- mean(df$NumberSelectedVars[df$MethodxNumberVars == i])
      dfAUC <- rbind(dfAUC, temp)
    }
    lims <- range(dfAUC$AUC,0.5) #set limits for AUC plots
    dfAUC$Codata <- "all"
    
    save(Res,df,dfAUC,file=fname)
  }
  
  #Fit model 4: all three co-data matrices, with positivity constraints----
  fname <- paste("Res4G",g,".Rdata",sep="")
  print(fname)
  
  if(run_withConstraints2){
    tic <- proc.time()[[3]]
    Res<-ecpc(Y=RespTrain, X=TrainData, #training data
              Z=Z_all[[g]], 
              paraPen = list(Z2=list(S1=S1_cor), Z3=list(S1=S1_pvals)),
              paraCon = list(Z2=Con_cor2, Z3=Con_pvals2),
              Y2=RespValidationNum, X2=ValidationData, #test data
              maxsel=maxSel) #maximum number of selected covariates for posterior selection
    toc <- proc.time()[[3]]; 
    print(toc - tic) #elapsed time
    Res$time <- toc-tic
    
    plot(Res, show = "priorweights", Z=Z_all[[g]], values=values[[g]])
    
    #concatenate results on predictions in data frame
    ntest <- length(RespValidation) #number of test samples in the validation data
    p<-dim(CoDataTrain)[1] #number of genes
    df<-data.frame("Ypred"= c(Res$Ypred, Res$Ypredridge, c(Res$YpredPost)),
                   "Truth" = rep(RespValidationNum, 2+length(maxSel)),
                   "Sample" = rep(1:ntest,2+length(maxSel)),
                   "NumberSelectedVars" = rep(c(sum(Res$beta!=0),p,maxSel),each=ntest),
                   "Method" = as.factor(rep(c("ecpc","ordinary.ridge",rep("ecpc",length(maxSel))),each=ntest)),
                   "Sparsity" = as.factor(rep(c("dense","dense", rep("sparse",length(maxSel))),each=ntest)),
                   "MethodxNumberVars" = as.factor(rep(c("ecpc","ordinary.ridge",
                                                         paste("ecpc",maxSel,"vars",sep="")),each=ntest)))
    
    #Compute AUC for each method and store in a data frame
    dfAUC <- data.frame()
    for (i in levels(df$MethodxNumberVars)) {
      temp <- data.frame("Method" = unique(df$Method[df$MethodxNumberVars==i]),
                         "Sparsity"=unique(df$Sparsity[df$MethodxNumberVars==i]),
                         "MethodxNumberVars"=i)
      rocpROC <- pROC::roc(predictor = df$Ypred[df$MethodxNumberVars == i], 
                           response = df$Truth[df$MethodxNumberVars == i], 
                           smooth = F, auc = T, levels = c(0, 1), direction = "<")
      temp$AUC <- rocpROC$auc[1]
      temp$NumberSelectedVars <- mean(df$NumberSelectedVars[df$MethodxNumberVars == i])
      dfAUC <- rbind(dfAUC, temp)
    }
    lims <- range(dfAUC$AUC,0.5) #set limits for AUC plots
    dfAUC$Codata <- "all"
    
    save(Res,df,dfAUC,file=fname)
  }
}

#Plots: general parameters----
wdth<-600
hght<-wdth*5/8
wdthpdf <- wdth/75
hghtpdf <- hght/75
ts <- 16 #basis text size in figures
ls <- 1.5 #basis line size in figures
ps <- 2 #basis point size in figures
sz <- 2 #point size
strk <- 1.5 #stroke size
palette <- "Dark2"

#Load data for plots----

#analyses with other co-data model settings
setting_names <- c("Positive+monotone constraints","Regular scale",
                   "Transformed scale","Positive constraints")
dfAUCall <- data.frame()
dfVk <- list(data.frame(),data.frame(),data.frame())
for(i in 1:length(setting_names)){
  for(g in 1:length(G_all)){
    fname <- paste("Res",i,"G",g,".Rdata",sep="")
    load(fname)
    dfAUC$Setting <- setting_names[i]
    dfAUC$Groups <- G_all[g]
    dfAUCall <- rbind(dfAUCall, dfAUC)
    
    indCodata <- attributes(Res$gamma)$codataSource
    #prior variance before truncation for each co-data source
    for(c in 1:length(Z_all[[g]])){
      vk <- as.vector(Res$tauglobal*(Z_all[[g]][[c]]%*%Res$gamma[indCodata==c]))
      if(c==3) vk <- as.vector(Res$tauglobal*(Z_all_rank[[g]][[c]]%*%Res$gamma[indCodata==c]))
      temp <- data.frame("Variable"=1:p,
                         "Priorvar"=vk,
                         "Priorvar*weight"=vk*Res$w[c])
      temp$Groups <- G_all[g]
      temp$Setting <- setting_names[i]
      temp$Gamma0 <- Res$gamma0
      temp$Tau2 <- Res$tauglobal
      temp$Sigma2 <- Res$sigmahat
      temp$Codata <- names(Z_all[[g]])[c]
      
      if(c==1){
        temp$values <- factor(CoDataTrain$RoepmanGenes,levels=c(0,1),
                              labels=c("out","in"))
      }else{
        temp$values <- values[[g]][[c]]
      }
      if(i%in%2:3){
        temp$Constraints <- F
      }else{
        temp$Constraints <- T
      }
      dfVk[[c]] <- rbind(dfVk[[c]],temp)
    }
  }
}
for(c in 1:3){
  dfVk[[c]]$Groups <- as.factor(dfVk[[c]]$Groups)
}
dfAUCall$Groups <- factor(dfAUCall$Groups)
dfAUCall$Setting <- factor(dfAUCall$Setting, levels=unique(dfAUCall$Setting)[c(2,3,4,1)],
                           labels=unique(dfAUCall$Setting)[c(2,3,4,1)])
dfAUCall$SettingNum <- factor(dfAUCall$Setting, levels=unique(dfAUCall$Setting)[c(2,3,4,1)],
                              labels=paste("Setting",1:4))
#Plot example----

plot(Res_example, show = "priorweights", Z=Z_example, values=values_example)
ggsave(filename="Plot_example.pdf",width= wdthpdf*1.5, height = hghtpdf*0.75)

#Plot contribution of each co-data source----

p1<-list()
for(c in 2:3){
  p1[[c-1]] <- ggplot(dfVk[[c]])+
    aes(x=values,y=Priorvar)+
    geom_line(aes(linetype=Groups),alpha=0.6,size=ls)+
    facet_grid(Codata~Setting)+
    labs(y="Prior variance")+
    theme_bw()+
    theme(axis.text.x=element_text(size=ts),
          axis.text.y=element_text(size=ts),
          axis.title.x=element_text(size=ts+2),
          axis.title.y=element_text(size=ts+2),
          legend.text=element_text(size=ts),
          legend.title=element_text(size=ts+2),
          legend.key.width = unit(2,"cm"),
          strip.text=element_text(size=ts))#,
  p1[[c-1]]
}
ggarrange(plotlist=p1, common.legend=TRUE, nrow=2, legend="right")

dfVk2 <- rbind(dfVk[[2]],dfVk[[3]])
dfVk2$SettingNum <- factor(dfVk2$Setting, levels=unique(dfVk2$Setting)[c(2,3,4,1)],
                           labels=paste("Setting",1:4))
figname <- paste("Estimates_LNM.pdf",sep="")
pdf(width = wdthpdf*1.2, height = hghtpdf,
    file = figname)
ggplot(dfVk2[dfVk2$Groups=="20",])+
  aes(x=values,y=Priorvar)+
  geom_line(aes(col=Setting),alpha=0.6,size=ls)+
  facet_grid(Constraints~Codata, scales="free", labeller="label_both")+
  labs(y="Prior variance", x="Continuous co-data variable")+
  guides(linetype=guide_legend(title="# Splines"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        legend.key.width = unit(2,"cm"),
        strip.text=element_text(size=ts))#,
dev.off()

figname <- paste("Estimates_LNM_groups.pdf",sep="")
pdf(width = wdthpdf*1.8, height = hghtpdf,
    file = figname)
ggplot(dfVk2)+
  aes(x=values,y=Priorvar)+
  geom_line(aes(linetype=Groups),alpha=0.6,size=ls)+
  facet_wrap(Codata~SettingNum, scales="free", nrow=2)+
  #facet_wrap(Codata~SettingNum, scales="free_x", nrow=2)+
  labs(y="Prior variance", x="Continuous co-data variable")+
  guides(linetype=guide_legend(title="# Splines"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        legend.key.width = unit(2,"cm"),
        strip.text=element_text(size=ts))#,
dev.off()

#Plot AUC----

lims <- range(dfAUCall$AUC,0.5) #set limits for AUC plots

figname <- paste("Prediction_LNM.pdf",sep="")
pdf(width = wdthpdf*1.8, height = hghtpdf,
    file = figname)
ggplot(dfAUCall[dfAUCall$Sparsity=="sparse",])+
  aes(x=NumberSelectedVars,y=AUC,col=Method)+
  geom_point(size=ps)+
  geom_line(size=ls,alpha=0.2)+
  facet_grid(Groups~SettingNum)+
  geom_hline(data=dfAUCall[dfAUCall$Sparsity=="dense",],
             aes(yintercept=AUC,col=Method), size=ls)+
  coord_cartesian(ylim=lims)+
  labs(x="# variables")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        #legend.key.width = unit(2,"cm"),
        legend.position = "bottom",
        strip.text=element_text(size=ts))#,
dev.off()







###############################################################################

#Session info----
sessionInfo()
