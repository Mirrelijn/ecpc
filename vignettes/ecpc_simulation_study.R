## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo=TRUE
)

## ---- eval=FALSE--------------------------------------------------------------
#  if(!requireNamespace("ggplot2")) install.packages("ggplot2")
#  if(!requireNamespace("ggpubr")) install.packages("ggpubr")
#  if(!requireNamespace("ecpc")) install.packages("ecpc")
#  if(!requireNamespace("squeezy")) install.packages("squeezy")
#  if(!requireNamespace("dplyr")) install.packages("dplyr")
#  if(!requireNamespace("RColorBrewer")) install.packages("RColorBrewer")
#  if(!requireNamespace("foreach")) install.packages("foreach")
#  if(!requireNamespace("doParallel")) install.packages("doParallel")
#  if(!requireNamespace("glmnet")) install.packages("glmnet")
#  if(!requireNamespace("mvtnorm")) install.packages("mvtnorm")
#  if(!requireNamespace("pROC")) install.packages("pROC")
#  if(!requireNamespace("devtools")) install.packages("devtools")
#  library(devtools)
#  if(!requireNamespace("fwelnet")) install_github("kjytay/fwelnet")
#  if(!requireNamespace("CoRF")) install_github("DennisBeest/CoRF")
#  if(!requireNamespace("ggh4x")) install_github("teunbrand/ggh4x")

## -----------------------------------------------------------------------------
pathResults <- "./Results_sim_study/"

## -----------------------------------------------------------------------------
runLinear <- FALSE
runGAM <- FALSE
runSCAM <- FALSE
runSCAMpmi <- FALSE
runAD <- FALSE

generateData <- FALSE #set to true if not generated before
nSim <- 50 #set to lower number for quicker run
runParallel <- FALSE

## -----------------------------------------------------------------------------
library(ecpc)
library(dplyr) #for data wrangling results
library(ggplot2) #for plotting results
library(RColorBrewer) #for plotting results
library(foreach) #for parallel computing
library(doParallel) #for parallel computing
library(scales)


if(runParallel){ #set to 1 to setup parallel backend to use many processors
  cores=detectCores()
  if(!("cl"%in%ls())){
    cl <- makeCluster(cores-1) #not to overload your computer
    registerDoParallel(cl)
  }
}

## -----------------------------------------------------------------------------
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
  save(Dat,AllData,file=paste(pathResults,"SimData.Rdata",sep=''))
}else{
  load(paste(pathResults,"SimData.Rdata",sep=''))
}


## -----------------------------------------------------------------------------

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
#plot(Z.all[[1]],Dat$beta^2)

#setting 2: informative, unequally spaced, information at edge
Z.all[["size.edge"]] <- abs(Dat$beta)
ZI.all[["size.edge"]] <- cbind(rep(1,p), Z.all[["size.edge"]])
#for generalised additive co-data model:
for(g in 1:length(G)){
  Zs.all[[g]][["size.edge"]] <- createZforSplines(values=Z.all[["size.edge"]], G=G[g], bdeg=3) 
}
#plot(Z.all[[2]],Dat$beta^2)

#setting 3: non-informative, transformed for equally spacing
Z.all[["noninformative+transformed"]] <- order(order(Z.all[["noninformative"]],
                                                     decreasing=FALSE),decreasing=FALSE)
ZI.all[["noninformative+transformed"]] <- cbind(rep(1,p), Z.all[["noninformative+transformed"]])
#for generalised additive co-data model:
for(g in 1:length(G)){
  Zs.all[[g]][["noninformative+transformed"]] <- createZforSplines(values=Z.all[["noninformative+transformed"]], G=G[g], bdeg=3) 
}
#plot(Z.all[[3]],Dat$beta^2)

#setting 4: informative, transformed for equally spacing
Z.all[["size.edge+transformed"]] <- order(order(abs(Dat$beta),decreasing=FALSE),decreasing=FALSE)
ZI.all[["size.edge+transformed"]] <- cbind(rep(1,p), Z.all[["size.edge+transformed"]])
for(g in 1:length(G)){
  Zs.all[[g]][["size.edge+transformed"]] <- createZforSplines(values=Z.all[["size.edge+transformed"]], G=G[g], bdeg=3) 
}
#plot(Z.all[[4]],Dat$beta^2)

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


## ---- eval=FALSE, class.source = 'fold-hide'----------------------------------
#  
#  #Fit ecpc on simulated data sets: linear co-data model----
#  fname <- paste(pathResults,"SimResAppNoteLinear",".Rdata",sep="") #(with intercept in ZI.all)
#  
#  if(runLinear){
#    #df <- data.frame()
#    #for(sim in 1:nSim){
#    finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
#                           .packages = c("ecpc")) %dopar% {
#                             df <- data.frame()
#                             for(setting in 1:4){
#                               for(method in "none"){#c("ML","fREML","GCV.Cp")){
#                                 tic<-proc.time()[[3]]
#                                 fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
#                                             Z=ZI.all[setting],
#                                             #bam.method=method,intrcpt.bam = F,
#                                             model="linear",maxsel=c(5,10,15,20),
#                                             Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
#                                             est_beta_method = "multiridge")
#                                 toc <- proc.time()[[3]]-tic
#  
#                                 vk <- (ZI.all[[setting]]%*%fit$gamma )*fit$tauglobal
#                                 #vk <- (matrix(Z.all[[setting]],p,1)%*%fit$gamma )*fit$tauglobal #without intercept
#                                 #plot(Z.all[[setting]],vk)
#  
#  
#                                 temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
#                                 temp$Z <- Z.all[[setting]]
#                                 temp$Time <- toc
#                                 temp$Covariate <- 1:p
#                                 temp$G <- 2
#                                 temp$setting <- setting
#                                 temp$bam.method <- method
#                                 temp$method <- "linear"
#                                 temp$MSEridge <- fit$MSEridge
#                                 temp$MSEecpc <- fit$MSEecpc
#                                 temp$Sim <- sim
#                                 temp$Codata <- names(Z.all)[setting]
#                                 temp$Transform <- F; if(setting>=3) temp$Transform <- T
#                                 temp$ZUntransformed <- Z.all[[setting]]
#                                 if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
#                                 df <- rbind(df,temp)
#                               }
#                             }
#                             list("df"=df)
#                           }
#    #str(finalMatrix)
#    df <- finalMatrix[[1]]; for(i in 2:nSim) df <- rbind(df,finalMatrix[[i]])
#    save(finalMatrix,df,file=fname)
#  }
#  

## ---- eval=FALSE, class.source = 'fold-hide'----------------------------------
#  #Fit ecpc on simulated data sets: generalised additive co-data model----
#  fname <- paste(pathResults,"SimResAppNoteGAM.Rdata",sep="")
#  print(fname)
#  
#  if(runGAM){
#    #df <- data.frame()
#    #for(sim in 1:nSim){
#    finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
#                           .packages = c("ecpc")) %dopar% {
#                             df <- data.frame()
#                             for(setting in 1:4){
#                               for(g in 1:length(G)){
#                                 for(method in c("splits","ML","fREML","GCV.Cp")){
#                                   tic<-proc.time()[[3]]
#  
#                                   if(method=="splits"){
#                                     fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
#                                                 Z=Zs.all[[g]][setting],
#                                                 paraPen=list(Z1=list(S1=S1[[g]])),
#                                                 intrcpt.bam = F,
#                                                 model="linear",maxsel=c(5,10,15,20),
#                                                 Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
#                                                 hypershrinkage = "ridge",
#                                                 est_beta_method = "multiridge")
#                                   }else{
#                                     fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
#                                                 Z=Zs.all[[g]][setting],
#                                                 paraPen=list(Z1=list(S1=S1[[g]])),
#                                                 bam.method=method,intrcpt.bam = F,
#                                                 model="linear",maxsel=c(5,10,15,20),
#                                                 Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
#                                                 est_beta_method = "multiridge")
#                                   }
#  
#                                   toc <- proc.time()[[3]]-tic
#  
#                                   vk <- (Zs.all[[g]][[setting]]%*%fit$gamma )*fit$tauglobal
#                                   #plot(Z.all[[setting]],vk)
#  
#  
#  
#                                   temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
#                                   temp$Z <- Z.all[[setting]]
#                                   temp$Time <- toc
#                                   temp$Covariate <- 1:p
#                                   temp$G <- G[g]
#                                   temp$setting <- setting
#                                   temp$bam.method <- method
#                                   temp$method <- "gam"
#                                   temp$MSEridge <- fit$MSEridge
#                                   temp$MSEecpc <- fit$MSEecpc
#                                   temp$Sim <- sim
#                                   temp$Codata <- names(Z.all)[setting]
#                                   temp$Transform <- F; if(setting>=3) temp$Transform <- T
#                                   temp$ZUntransformed <- Z.all[[setting]]
#                                   if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
#                                   df <- rbind(df,temp)
#                                 }
#                               }
#                             }
#                             list("df"=df)
#                           }
#    #str(finalMatrix)
#    df <- finalMatrix[[1]]; for(i in 2:nSim) df <- rbind(df,finalMatrix[[i]])
#    save(finalMatrix,df,file=fname)
#  }
#  

## ---- eval=FALSE, class.source = 'fold-hide'----------------------------------
#  #Fit ecpc on simulated data sets: positive constrained generalised additive co-data model----
#  fname <- paste(pathResults,"SimResAppNoteSCAMp.Rdata",sep="")
#  print(fname)
#  
#  if(runSCAM){
#    #df <- data.frame()
#    #for(sim in 1:nSim){
#    finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
#                           .packages = c("ecpc")) %dopar% {
#                             df <- data.frame()
#                             for(setting in 1:4){
#                               for(g in 1:length(G)){
#                                 for(method in c("splits")){
#                                   tic<-proc.time()[[3]]
#  
#                                   fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
#                                               Z=Zs.all[[g]][setting],
#                                               paraPen=list(Z1=list(S1=S1[[g]])),
#                                               paraCon = list(Z1=Con.p[[g]]),
#                                               intrcpt.bam = F,
#                                               model="linear",maxsel=c(5,10,15,20),
#                                               Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
#                                               est_beta_method = "multiridge")
#                                   toc <- proc.time()[[3]]-tic
#  
#                                   vk <- (Zs.all[[g]][[setting]]%*%fit$gamma )*fit$tauglobal
#                                   plot(Z.all[[setting]],vk)
#  
#  
#  
#                                   temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
#                                   temp$Z <- Z.all[[setting]]
#                                   temp$Time <- toc
#                                   temp$Covariate <- 1:p
#                                   temp$G <- G[g]
#                                   temp$setting <- setting
#                                   temp$bam.method <- method
#                                   temp$method <- "scam.p"
#                                   temp$MSEridge <- fit$MSEridge
#                                   temp$MSEecpc <- fit$MSEecpc
#                                   temp$Sim <- sim
#                                   temp$Codata <- names(Z.all)[setting]
#                                   temp$Transform <- F; if(setting>=3) temp$Transform <- T
#                                   temp$ZUntransformed <- Z.all[[setting]]
#                                   if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
#                                   df <- rbind(df,temp)
#                                 }
#                               }
#                             }
#                             list("df"=df)
#                           }
#    #str(finalMatrix)
#    df <- finalMatrix[[1]]; for(i in 2:nSim) df <- rbind(df,finalMatrix[[i]])
#    save(finalMatrix,df,file=fname)
#  }
#  

## ---- eval=FALSE, class.source = 'fold-hide'----------------------------------
#  #Fit ecpc on simulated data sets: positive+monotone increasing constrained generalised additive co-data model----
#  fname <- paste(pathResults,"SimResAppNoteSCAMpmi.Rdata",sep="")
#  print(fname)
#  
#  if(runSCAMpmi){
#    #df <- data.frame()
#    #for(sim in 1:nSim){
#    finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
#                           .packages = c("ecpc")) %dopar% {
#                             df <- data.frame()
#                             for(setting in 1:4){
#                               for(g in 1:length(G)){
#                                 for(method in c("splits")){
#                                   tic<-proc.time()[[3]]
#  
#                                   fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
#                                               Z=Zs.all[[g]][setting],
#                                               paraPen=list(Z1=list(S1=S1[[g]])),
#                                               paraCon = list(Z1=Con.pmi[[g]]),
#                                               intrcpt.bam = F,
#                                               model="linear",maxsel=c(5,10,15,20),
#                                               Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
#                                               est_beta_method = "multiridge")
#                                   toc <- proc.time()[[3]]-tic
#  
#                                   vk <- (Zs.all[[g]][[setting]]%*%fit$gamma )*fit$tauglobal
#                                   #plot(Z.all[[setting]],vk)
#  
#  
#  
#                                   temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
#                                   temp$Z <- Z.all[[setting]]
#                                   temp$Time <- toc
#                                   temp$Covariate <- 1:p
#                                   temp$G <- G[g]
#                                   temp$setting <- setting
#                                   temp$bam.method <- method
#                                   temp$method <- "scam.pmi"
#                                   temp$MSEridge <- fit$MSEridge
#                                   temp$MSEecpc <- fit$MSEecpc
#                                   temp$Sim <- sim
#                                   temp$Codata <- names(Z.all)[setting]
#                                   temp$Transform <- F; if(setting>=3) temp$Transform <- T
#                                   temp$ZUntransformed <- Z.all[[setting]]
#                                   if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
#                                   df <- rbind(df,temp)
#                                 }
#                               }
#                             }
#                             list("df"=df)
#                           }
#    #str(finalMatrix)
#    df <- finalMatrix[[1]]; for(i in 2:nSim) df <- rbind(df,finalMatrix[[i]])
#    save(finalMatrix,df,file=fname)
#  }

## ---- eval=FALSE, class.source = 'fold-hide'----------------------------------
#  #Fit ecpc on simulated data sets: adaptive discretisation co-data model----
#  fname <- paste(pathResults,"SimResAppNoteAD.Rdata",sep="")
#  print(fname)
#  
#  if(runAD){
#    #df <- data.frame()
#    #for(sim in 1:nSim){
#    finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
#                           .packages = c("ecpc")) %dopar% {
#                             df <- data.frame()
#                             for(setting in 1:2){
#                               #for(g in 1:length(G)){
#                               for(method in c("splits")){
#                                 tic<-proc.time()[[3]]
#  
#                                 fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
#                                             groupsets=groupsets[setting],
#                                             groupsets.grouplvl = hierarchy.grouplevel[setting],
#                                             hypershrinkage = "hierLasso,ridge",
#                                             model="linear",maxsel=c(5,10,15,20),
#                                             Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
#                                             est_beta_method = "multiridge")
#                                 toc <- proc.time()[[3]]-tic
#  
#                                 vk <- as.vector(Z.AD[[setting]]%*%fit$gamma*fit$tauglobal)
#                                 #plot(Z.all[[setting]],vk)
#  
#  
#                                 temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
#                                 temp$Z <- Z.all[[setting]]
#                                 temp$Time <- toc
#                                 temp$Covariate <- 1:p
#                                 temp$G <- length(groupsets[[setting]])
#                                 temp$setting <- setting
#                                 temp$bam.method <- method
#                                 temp$method <- "AD"
#                                 temp$MSEridge <- fit$MSEridge
#                                 temp$MSEecpc <- fit$MSEecpc
#                                 temp$Sim <- sim
#                                 temp$Codata <- names(Z.all)[setting]
#                                 temp$Transform <- F; if(setting>=3) temp$Transform <- T
#                                 temp$ZUntransformed <- Z.all[[setting]]
#                                 if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
#                                 df <- rbind(df,temp)
#                               }
#                               #}
#                             }
#                             list("df"=df)
#                           }
#    #str(finalMatrix)
#    df <- finalMatrix[[1]]; for(i in 2:nSim) df <- rbind(df,finalMatrix[[i]])
#    save(finalMatrix,df,file=fname)
#  }
#  

## -----------------------------------------------------------------------------
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
colpal <- "Dark2"
colsfill <- brewer.pal(3,"Dark2")[1:2]
colsAUC <- brewer.pal(10,"RdYlBu")
colsAUC <- seq_gradient_pal(low="black",high="white")((0:4)/5)



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


## ---- fig.width=12, fig.height=16---------------------------------------------
#Figure estimated variance vs continuous co-data for all simulations, fixed G and fixed bam.method----
g <- G[1]
ggplot(dfEst[dfEst$G%in%c(1,2,g,7) & dfEst$Transform==F &
               dfEst$bam.method%in%c("ML","none","splits"),])+# & dfEst$method=="ML",])+
  aes(x=Z)+
  #facet_wrap(Codata~method,scales="free_x",nrow = 2)+
  facet_grid(method~Codata,scales="free_x")+
  geom_point(aes(x=Z,y=truevk),alpha=0.3,col="grey20")+
  geom_ribbon(aes(ymin=q05Vk,ymax=q95Vk),linetype=0,alpha=0.1)+ 
  geom_ribbon(aes(ymin=q25Vk,ymax=q75Vk),linetype=0,alpha=0.2)+ 
  geom_line(aes(y=q50Vk),linewidth=1)+ #median
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


## ---- fig.width=12, fig.height=8----------------------------------------------
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

## ---- fig.width=12, fig.height=8----------------------------------------------
#Figure estimated variance vs continuous co-data for 1 run, different bam.methods----
g <- G[1]
Sim <- 1
ggplot(dfGAMs[dfGAMs$Sim==Sim & dfGAMs$Transform==F & dfGAMs$G%in%c(1,G) ,])+
  aes(x=Z,col=bam.method)+
  facet_grid(.~Codata,scales="free_x")+
  geom_point(aes(x=Z,y=truevk),alpha=0.1,col="black")+
  geom_line(aes(y=vkfit,linetype=G),linewidth=1.5)+
  labs(y="Prior variance",x="Continuous co-data variable")+
  scale_color_manual(values=colsAUC[1:4])+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,


## ---- fig.width=12, fig.height=12---------------------------------------------
#Figure estimated variance vs continuous co-data for all simulations, gam only, fixed G and different bam.method----
g <- G[1]
ggplot(dfEstGAMs[dfEstGAMs$G%in%c(1,2,g) & dfEstGAMs$Transform==F ,])+
  aes(x=Z)+
  facet_grid(bam.method~Codata,scales="free_x")+
  geom_point(aes(x=Z,y=truevk),alpha=0.3,col="black")+
  geom_ribbon(aes(ymin=q05Vk,ymax=q95Vk),linetype=0,alpha=0.1)+ 
  geom_ribbon(aes(ymin=q25Vk,ymax=q75Vk),linetype=0,alpha=0.2)+ 
  geom_line(aes(y=q50Vk),linewidth=1)+ #median
  labs(y="Prior variance",x="Continuous co-data variable")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,

## ---- fig.width=12, fig.height=8----------------------------------------------
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

ggplot(temp2[temp2$Transform==F,])+#& !(dfPred$Codata=="noninformative"),])+
  aes(x=bam.method,y=MSE)+
  #geom_boxplot(data=temp2[temp2$method=="ridge",],aes(fill=method))+ #ridge prediction performance
  geom_boxplot(aes(fill=G))+
  facet_grid(.~Codata2)+
  #scale_y_log10()+-
  labs(y="MSE",x="")+
  scale_fill_manual(values=rev(rev(colsAUC)[1:3]))+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts,angle=30,vjust=1,hjust=1),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,

## -----------------------------------------------------------------------------
#Table average computing times + standard deviation----
timeTable <- dfAll[dfAll$Covariate==1,] %>% group_by(method,bam.method,G) %>% 
  summarise(meanTime=mean(Time),sdTime=sd(Time)) %>% ungroup()
print(timeTable)
#library("writexl")
#fname <- paste(pathFigures,"TimeTableRaw.xlsx")
#write_xlsx(timeTable,fname)

## -----------------------------------------------------------------------------
run_ecpc <- FALSE
run_glmnet <- FALSE
run_fwelnet <- FALSE

generateData <- FALSE #set to true if not generated before
nSim <- 50 #set to lower number for quicker run
runParallel <- FALSE

## -----------------------------------------------------------------------------
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
q <- 5/6 #proportion zeros
maxsel <- c(rep(2:10),10*2:29)
alp_range <- seq(0,1,length.out = length(maxsel))

## -----------------------------------------------------------------------------
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
#hist(Dat$lp2)
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
  }
  save(Dat,AllData,file=paste(pathResults,"SimDataSparse",setting_q,rho,".Rdata",sep=""))
}else{
  load(paste(pathResults,"SimDataSparse",setting_q,rho,".Rdata",sep=""))
}


## -----------------------------------------------------------------------------
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
#plot(Z.all[[1]],Dat$beta^2)

#setting 2: informative, unequally spaced, information at edge
set.seed(101010)
Z.all[["size.edge"]] <- abs(Dat$beta)+rnorm(p,sd=sd(Dat$beta)/10)
ZI.all[["size.edge"]] <- cbind(rep(1,p), Z.all[["size.edge"]])
#for generalised additive co-data model:
for(g in 1:length(G)){
  Zs.all[[g]][["size.edge"]] <- createZforSplines(values=Z.all[["size.edge"]], G=G[g], bdeg=3) 
}
#plot(Z.all[[2]],Dat$beta^2)

#setting 3: informative, unequally spaced, information at 2 edges
set.seed(101010)
Z.all[["size.edge2"]] <- sign(Dat$beta)*abs(Dat$beta)+rnorm(p,sd=sd(Dat$beta)/10)
ZI.all[["size.edge2"]] <- cbind(rep(1,p), Z.all[["size.edge2"]])
#for generalised additive co-data model:
for(g in 1:length(G)){
  Zs.all[[g]][["size.edge2"]] <- createZforSplines(values=Z.all[["size.edge2"]], G=G[g], bdeg=3) 
}
#plot(Z.all[[3]],Dat$beta^2)

#fwelnet
Z_fwelnet <- list(matrix(Z.all[["noninformative"]],nrow = p),
                  matrix(Z.all[["size.edge"]],nrow = p),
                  matrix(Z.all[["size.edge2"]],nrow = p))

## ---- eval=FALSE, class.source = 'fold-hide'----------------------------------
#  #Fit ecpc+posterior selection on simulated data sets: generalised additive co-data model----
#  fname <- paste(pathResults,"SimResAppNoteSparse_GAM",setting_q,rho,".Rdata",sep="")
#  #fname <- paste(pathResults,"SimResAppNoteGAMSparse_splits",q,rho,".Rdata",sep="") #hypershrinkage=ridge
#  print(fname)
#  
#  if(run_ecpc){
#    #df <- data.frame()
#    #for(sim in 1:nSim){
#    finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
#                           .packages = c("ecpc", "squeezy")) %dopar% {
#                             df <- data.frame()
#                             dfSelect <- data.frame()
#                             for(setting in 1:3){
#                               for(g in 1){
#                                 for(method in c("ML")){
#                                   tic<-proc.time()[[3]]
#  
#                                   if(method=="splits"){
#                                     fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
#                                                 Z=Zs.all[[g]][setting],
#                                                 paraPen=list(Z1=list(S1=S1[[g]])),
#                                                 intrcpt.bam = F,
#                                                 model=model,maxsel=maxsel,
#                                                 Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
#                                                 hypershrinkage = "ridge",
#                                                 est_beta_method = "multiridge")
#                                   }else{
#                                     fit <- ecpc(Y=AllData[[sim]]$Y,X=AllData[[sim]]$Xctd,
#                                                 Z=Zs.all[[g]][setting],
#                                                 paraPen=list(Z1=list(S1=S1[[g]])),
#                                                 bam.method=method,intrcpt.bam = F,
#                                                 model=model,maxsel=maxsel,
#                                                 Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd,
#                                                 est_beta_method = "multiridge")
#                                   }
#  
#  
#                                   toc <- proc.time()[[3]]-tic
#  
#                                   vk <- (Zs.all[[g]][[setting]]%*%fit$gamma )*fit$tauglobal
#                                   #plot(Z.all[[setting]],vk)
#  
#                                   temp<-data.frame(vkfit=vk,truebeta=Dat$beta,truevk=Dat$beta^2)
#                                   temp$Z <- Z.all[[setting]]
#                                   temp$Time <- toc
#                                   temp$Covariate <- 1:p
#                                   temp$G <- G[g]
#                                   temp$setting <- setting
#                                   temp$bam.method <- method
#                                   temp$method <- "gam"
#                                   temp$MSEridge <- fit$MSEridge
#                                   temp$MSEecpc <- fit$MSEecpc
#                                   temp$Sim <- sim
#                                   temp$Codata <- names(Z.all)[setting]
#                                   temp$Transform <- F; if(setting>=3) temp$Transform <- T
#                                   temp$ZUntransformed <- Z.all[[setting]]
#                                   temp$q <- q
#                                   temp$rho <- rho
#                                   if(setting>=3) temp$ZUntransformed <- Z.all[[setting-2]]
#                                   df <- rbind(df,temp)
#  
#                                   #sensitivity = true positive rate = true positives/total positives
#                                   sensitivity <- apply(fit$betaPost,2,function(x){
#                                     ind_nonzero_est <- which(x!=0)
#                                     TPR <-sum(ind_nonzero_est%in%ind_nonzero)/length(ind_nonzero)
#                                     return(TPR)
#                                   })
#                                   #precision = true negative rate = true negatives/total negatives
#                                   precision <- apply(fit$betaPost,2,function(x){
#                                     ind_zero_est <- which(x==0)
#                                     TNR <-sum(ind_zero_est%in%ind_zero)/length(ind_zero)
#                                     return(TNR)
#                                   })
#                                   lp_train <- AllData[[sim]]$Xctd%*%fit$betaPost
#                                   MSEtrain <- apply(lp_train,2,function(lp) mean((lp-AllData[[sim]]$Y)^2))
#  
#                                   temp<-data.frame("Tuningparam"=maxsel, sensitivity, precision)
#                                   temp$TypeTuning <- "#params"
#                                   temp$Time <- toc
#                                   temp$G <- G[g]
#                                   temp$setting <- setting
#                                   temp$bam.method <- method
#                                   temp$method <- "gam"
#                                   temp$MSEridge <- fit$MSEridge
#                                   temp$MSEecpc <- c(fit$MSEPost)
#                                   temp$MSEtrain <- c(MSEtrain)
#                                   if(is.null(fit$MSEPost)) temp$MSEecpc <- rep(NaN,length(maxsel))
#                                   temp$Sim <- sim
#                                   temp$Codata <- names(Z.all)[setting]
#                                   temp$Transform <- F; if(setting>=3) temp$Transform <- T
#                                   temp$q <- q
#                                   temp$setting_q <- setting_q
#                                   temp$rho <- rho
#                                   dfSelect <- rbind(dfSelect, temp)
#  
#                                   #fit squeezy
#                                   #Use squeezy function to transform estimated ridge penalties to elastic net
#                                   #penalties
#                                   notinf <- fit$penalties!=Inf
#                                   fit.EN <- lapply(alp_range,function(alp){
#                                     fit.EN <- squeezy(Y=AllData[[sim]]$Y, X=AllData[[sim]]$Xctd[,notinf],
#                                                       groupset=lapply(1:sum(notinf), function(x) x),
#                                                       alpha=alp,
#                                                       Y2=AllData[[sim]]$Y2,X2=AllData[[sim]]$X2ctd[,notinf],
#                                                       lambdas=fit$penalties[notinf],
#                                                       sigmasq=fit$sigmahat)})
#                                   beta_squeezy <- sapply(1:length(alp_range),function(i){
#                                     betas <- rep(0,p)
#                                     betas[notinf] <- fit.EN[[i]]$betaApprox
#                                     return(betas)
#                                   })
#                                   MSE_squeezy <- sapply(1:length(alp_range),function(i){
#                                     fit.EN[[i]]$MSEApprox
#                                   })
#                                   lp_train <- AllData[[sim]]$Xctd%*%beta_squeezy
#                                   MSEtrain <- apply(lp_train,2,function(lp) mean((lp-AllData[[sim]]$Y)^2))
#  
#                                   #sensitivity = true positive rate = true positives/total positives
#                                   sensitivity <- apply(beta_squeezy,2,function(x){
#                                     ind_nonzero_est <- which(x!=0)
#                                     TPR <-sum(ind_nonzero_est%in%ind_nonzero)/length(ind_nonzero)
#                                     return(TPR)
#                                   })
#                                   #precision = true negative rate = true negatives/total negatives
#                                   precision <- apply(beta_squeezy,2,function(x){
#                                     ind_zero_est <- which(x==0)
#                                     TNR <-sum(ind_zero_est%in%ind_zero)/length(ind_zero)
#                                     return(TNR)
#                                   })
#  
#                                   temp<-data.frame("Tuningparam"=alp_range, sensitivity, precision)
#                                   temp$TypeTuning <- "alpha"
#                                   temp$Time <- toc
#                                   temp$G <- G[g]
#                                   temp$setting <- setting
#                                   temp$bam.method <- method
#                                   temp$method <- "gam"
#                                   temp$MSEridge <- fit$MSEridge
#                                   temp$MSEecpc <- MSE_squeezy
#                                   temp$MSEtrain <- c(MSEtrain)
#                                   temp$Sim <- sim
#                                   temp$Codata <- names(Z.all)[setting]
#                                   temp$Transform <- F; if(setting>=3) temp$Transform <- T
#                                   temp$q <- q
#                                   temp$setting_q <- setting_q
#                                   temp$rho <- rho
#                                   dfSelect <- rbind(dfSelect,temp)
#  
#                                 }
#                               }
#                             }
#                             list("df"=df,"dfSelect"=dfSelect)
#                           }
#    #str(finalMatrix)
#    df2 <- lapply(1:nSim,function(i) finalMatrix[i,1][[1]])
#    dfSelect2 <- lapply(1:nSim,function(i) finalMatrix[i,2][[1]])
#    df <- df2[[1]]; for(i in 2:nSim) df <- rbind(df,df2[[i]])
#    dfSelect <- dfSelect2[[1]]; for(i in 2:nSim) dfSelect <- rbind(dfSelect,dfSelect2[[i]])
#    save(finalMatrix,df,dfSelect,file=fname)
#  }

## ---- eval=FALSE, class.source = 'fold-hide'----------------------------------
#  #Fit glmnet on simulated data sets----
#  fname <- paste(pathResults,"SimResAppNoteSparse_glmnet",setting_q,rho,".Rdata",sep="")
#  #fname <- paste(pathResults,"SimResAppNoteGAMSparse_splits",q,rho,".Rdata",sep="") #hypershrinkage=ridge
#  print(fname)
#  
#  if(run_glmnet){
#    #df <- data.frame()
#    #for(sim in 1:nSim){
#    finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
#                           .packages = c("glmnet")) %dopar% {
#                             dfSelect <- data.frame()
#                             for(setting in 1:3){
#                               for(alp in alp_range){
#                                 tic<-proc.time()[[3]]
#  
#                                 fit.glmnet <- glmnet::cv.glmnet(y=AllData[[sim]]$Y,x=AllData[[sim]]$Xctd,
#                                                                 family=fam,alpha=alp)
#                                 beta.glmnet <- coef(fit.glmnet,s="lambda.min", exact=TRUE)
#                                 Ypred.glmnet <- c(predict(fit.glmnet, newx = AllData[[sim]]$X2ctd,
#                                                           s = "lambda.min", type="response", exact=TRUE))
#                                 MSE.glmnet <- mean((Ypred.glmnet-AllData[[sim]]$Y2)^2)
#                                 lp_train <- AllData[[sim]]$Xctd%*%beta.glmnet[-1] + beta.glmnet[1]
#                                 MSEtrain <- mean((lp_train-AllData[[sim]]$Y)^2)
#                                 toc <- proc.time()[[3]]-tic
#  
#                                 #sensitivity = true positive rate = true positives/total positives
#                                 ind_nonzero_est <- which(beta.glmnet[-1]!=0)
#                                 sensitivity <- sum(ind_nonzero_est%in%ind_nonzero)/length(ind_nonzero)
#                                 #precision = true negative rate = true negatives/total negatives
#                                 ind_zero_est <- which(beta.glmnet[-1]==0)
#                                 precision <- sum(ind_zero_est%in%ind_zero)/length(ind_zero)
#  
#                                 temp<-data.frame("Tuningparam"=alp, sensitivity, precision)
#                                 temp$TypeTuning <- "alpha"
#                                 temp$Time <- toc
#                                 temp$G <- 1
#                                 temp$setting <- setting
#                                 temp$bam.method <- "none"
#                                 temp$method <- "glmnet"
#                                 temp$MSEridge <- NaN
#                                 temp$MSEecpc <- MSE.glmnet
#                                 temp$MSEtrain <- MSEtrain
#                                 temp$Sim <- sim
#                                 temp$Codata <- names(Z.all)[setting]
#                                 temp$Transform <- F; if(setting>=3) temp$Transform <- T
#                                 temp$q <- q
#                                 temp$setting_q <- setting_q
#                                 temp$rho <- rho
#                                 dfSelect <- rbind(dfSelect,temp)
#  
#                               }
#                             }
#                             list("dfSelect"=dfSelect)
#                           }
#    #str(finalMatrix)
#    dfSelect <- finalMatrix[[1]]; for(i in 2:nSim) dfSelect <- rbind(dfSelect,finalMatrix[[i]])
#    save(finalMatrix,dfSelect,file=fname)
#  }

## ---- eval=FALSE, class.source = 'fold-hide'----------------------------------
#  #Fit fwelnet on simulated data sets----
#  fname <- paste(pathResults,"SimResAppNoteSparse_fwelnet",setting_q,rho,".Rdata",sep="")
#  print(fname)
#  
#  if(run_fwelnet){
#    #df <- data.frame()
#    #for(sim in 1:nSim){
#    finalMatrix <- foreach(sim=1:nSim, .combine=rbind,
#                           .packages = c("fwelnet")) %dopar% {
#                             dfSelect <- data.frame()
#                             for(setting in 1:3){
#                               for(alp in alp_range){
#                                 tic<-proc.time()[[3]]
#  
#                                 res.fwEN <-  cv.fwelnet(y=AllData[[sim]]$Y,x=AllData[[sim]]$Xctd,
#                                                         z=Z_fwelnet[[setting]],
#                                                         alpha=alp,family=fam,standardize=F)
#  
#                                 ind.minlam <- which(res.fwEN$lambda==res.fwEN$lambda.min)
#  
#                                 betafwEN <- res.fwEN$glmfit$beta[,ind.minlam]
#                                 a0fwEN <- res.fwEN$glmfit$a0[ind.minlam]
#                                 Ypred_fwelnet <- AllData[[sim]]$X2ctd%*%betafwEN+a0fwEN
#                                 MSE_fwelnet <- mean((Ypred_fwelnet-AllData[[sim]]$Y2)^2)
#                                 lp_train <- AllData[[sim]]$Xctd%*%betafwEN + a0fwEN
#                                 MSEtrain <- mean((lp_train-AllData[[sim]]$Y)^2)
#  
#                                 toc <- proc.time()[[3]]-tic
#  
#                                 #sensitivity = true positive rate = true positives/total positives
#                                 ind_nonzero_est <- which(betafwEN!=0)
#                                 sensitivity <- sum(ind_nonzero_est%in%ind_nonzero)/length(ind_nonzero)
#                                 #precision = true negative rate = true negatives/total negatives
#                                 ind_zero_est <- which(betafwEN==0)
#                                 precision <- sum(ind_zero_est%in%ind_zero)/length(ind_zero)
#  
#                                 temp<-data.frame("Tuningparam"=alp, sensitivity, precision)
#                                 temp$TypeTuning <- "alpha"
#                                 temp$Time <- toc
#                                 temp$G <- 1
#                                 temp$setting <- setting
#                                 temp$bam.method <- "none"
#                                 temp$method <- "fwelnet"
#                                 temp$MSEridge <- NaN
#                                 temp$MSEecpc <- MSE_fwelnet
#                                 temp$MSEtrain <- c(MSEtrain)
#                                 temp$Sim <- sim
#                                 temp$Codata <- names(Z.all)[setting]
#                                 temp$Transform <- F; if(setting>=3) temp$Transform <- T
#                                 temp$q <- q
#                                 temp$setting_q <- setting_q
#                                 temp$rho <- rho
#                                 dfSelect <- rbind(dfSelect,temp)
#  
#                               }
#                             }
#                             list("dfSelect"=dfSelect)
#                           }
#    #str(finalMatrix)
#    dfSelect <- finalMatrix[[1]]; for(i in 2:nSim) dfSelect <- rbind(dfSelect,finalMatrix[[i]])
#    save(finalMatrix,dfSelect,file=fname)
#  }

## -----------------------------------------------------------------------------
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
  summarise(meanMSE=mean(MSEecpc,na.rm=TRUE),method2="EN",
            q50MSE = quantile(MSEecpc,0.5,na.rm=TRUE),
            q95MSE = quantile(MSEecpc,0.95,na.rm=TRUE),q05MSE = quantile(MSEecpc,0.05,na.rm=TRUE),
            q75MSE = quantile(MSEecpc,0.75,na.rm=TRUE),q25MSE = quantile(MSEecpc,0.25,na.rm=TRUE),
            meanMSEtrain=mean(MSEtrain,na.rm=TRUE),method2="EN",
            q50MSEtrain = quantile(MSEtrain,0.5,na.rm=TRUE),
            q95MSEtrain = quantile(MSEtrain,0.95,na.rm=TRUE),q05MSEtrain = quantile(MSEtrain,0.05,na.rm=TRUE),
            q75MSEtrain = quantile(MSEtrain,0.75,na.rm=TRUE),q25MSEtrain = quantile(MSEtrain,0.25,na.rm=TRUE),
            meansensitivity=mean(sensitivity,na.rm=TRUE),
            q50sens = quantile(sensitivity,0.5,na.rm=T),
            q95sens = quantile(sensitivity,0.95,na.rm=T),q05sens = quantile(sensitivity,0.05,na.rm=T),
            q75sens = quantile(sensitivity,0.75,na.rm=T),q25sens = quantile(sensitivity,0.25,na.rm=T),
            meanprecision=mean(precision,na.rm=TRUE),
            q50prec = quantile(precision,0.5,na.rm=T),
            q95prec = quantile(precision,0.95,na.rm=T),q05prec = quantile(precision,0.05,na.rm=T),
            q75prec = quantile(precision,0.75,na.rm=T),q25prec = quantile(precision,0.25,na.rm=T)) %>% ungroup()


## ---- fig.width=12, fig.height=6----------------------------------------------
#Figure MSE for all simulations, tuning parameter 2----
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


## ---- fig.width=12, fig.height=6----------------------------------------------
#Figure sensitivity vs precision for all simulations----
ggplot(dfPred)+# & dfEst$method=="ML",])+
  aes(col=method, shape=TypeTuning2)+
  facet_grid(.~Codata)+
  geom_line(aes(y=meansensitivity,x=1-meanprecision),linewidth=ls, alpha=0.2)+ #mean
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

## -----------------------------------------------------------------------------
run_n <- FALSE
run_p <- FALSE
generateData <- FALSE

## -----------------------------------------------------------------------------
#Load libraries-----
library(peakRAM) #for keeping track of time and peak memory
library(ecpc)
library(dplyr) #for data wrangling results
library(ggplot2) #for plotting results
library(ggpubr)
library(RColorBrewer) #for plotting results

## -----------------------------------------------------------------------------
# Simulate toy data ------
p_all <- c(1000, 2000, 5000, 10000, 20000, 50000, 100000)
n_all <- c(50, 100, 200, 500, 1000)

n2<-100 #sample size test data set
sigma <- 1
tauglobal <- 0.1 #prior variance

#simulate all betas i.i.d. from beta_k~N(mean=0,sd=sqrt(0.1)):
set.seed(5463)
if(generateData){
  AllData_p <- list()
  for(i in 1:length(p_all)){
    p<-p_all[i] #number of covariates
    n<-200 #sample size training data set
    AllData_p[[i]] <- list()
    AllData_p[[i]]$beta <- rnorm(p,sd=sqrt(tauglobal))
    AllData_p[[i]]$Xctd <- matrix(rnorm(n*p) ,n,p)
    AllData_p[[i]]$X2ctd <- matrix(rnorm(n*p) ,n,p)
    AllData_p[[i]]$Y <- rnorm(n,mean=c(AllData_p[[i]]$Xctd%*%AllData_p[[i]]$beta),sd=sigma)
    AllData_p[[i]]$Y2 <- rnorm(n2,mean=c(AllData_p[[i]]$X2ctd%*%AllData_p[[i]]$beta),sd=sigma)
  }
  AllData_n <- list()
  
  p<-5000 #number of covariates
  betas <- rnorm(p,sd=sqrt(tauglobal))
  for(i in 1:length(n_all)){
    n<-n_all[i] #sample size training data set
    AllData_n[[i]] <- list()
    AllData_n[[i]]$beta <- betas
    AllData_n[[i]]$Xctd <- matrix(rnorm(n*p) ,n,p)
    AllData_n[[i]]$X2ctd <- matrix(rnorm(n*p) ,n,p)
    AllData_n[[i]]$Y <- rnorm(n,mean=c(AllData_n[[i]]$Xctd%*%AllData_n[[i]]$beta),sd=sigma)
    AllData_n[[i]]$Y2 <- rnorm(n2,mean=c(AllData_n[[i]]$X2ctd%*%AllData_n[[i]]$beta),sd=sigma)
  }
  save(AllData_p,AllData_n,file=paste(pathResults,"SimData_p_n.Rdata",sep=''))
}else{
  load(paste(pathResults,"SimData_p_n.Rdata",sep=''))
}


## ---- eval=FALSE, class.source = 'fold-hide'----------------------------------
#  
#  #Run for p=5000, different n and check time and peak memory----
#  if(run_n){
#    fname <- paste(pathResults, 'res_n.Rdata', sep='')
#    df_n = data.frame("Method"=c(), "Time"=c(), "Peak_memory"=c(),
#                      "n" = c(), "p"=c())
#    p=5000
#    #Compute co-data for p=5000
#    Z <- abs(AllData_n[[1]]$beta)
#    ZI <- cbind(rep(1,p),abs(AllData_n[[1]]$beta))
#    Zs <- createZforSplines(values=Z, G=20, bdeg=3)
#    S1 <- createS(orderPen=2, G=20) #create 2nd order difference penalty matrix (same for all co-data)
#    Con.p <- createCon(G=20, shape="positive")
#  
#    for(i in 1:length(n_all)){
#      n = n_all[i]
#  
#      #ecpc linear co-data model
#      mem <- peakRAM({
#        fit <- ecpc(Y=AllData_n[[i]]$Y,X=AllData_n[[i]]$Xctd,
#                  Z=list(ZI),
#                  model="linear",postselection=FALSE,est_beta_method = "multiridge")
#      })
#      temp <- data.frame("Method"="ecpc_linear", "Time"=mem[[2]], "Peak_memory"=mem[[4]],
#                         "n"=dim(AllData_n[[i]]$Xctd)[1], p=dim(AllData_n[[i]]$Xctd)[2])
#      df_n <- rbind(df_n, temp)
#  
#      #ecpc GAM co-data model
#      mem <- peakRAM({
#        fit <- ecpc(Y=AllData_n[[i]]$Y,X=AllData_n[[i]]$Xctd,
#                    Z=list(Zs),
#                    paraPen=list(Z1=list(S1=S1)),
#                    bam.method="ML", intrcpt.bam = F,
#                    model="linear",postselection=FALSE,
#                    est_beta_method = "multiridge")
#      })
#      temp <- data.frame("Method"="ecpc_gam", "Time"=mem[[2]], "Peak_memory"=mem[[4]],
#                         "n"=dim(AllData_n[[i]]$Xctd)[1], p=dim(AllData_n[[i]]$Xctd)[2])
#      df_n <- rbind(df_n, temp)
#  
#      #ecpc SCAM co-data model
#      mem <- peakRAM({
#        fit <- ecpc(Y=AllData_n[[i]]$Y,X=AllData_n[[i]]$Xctd,
#                    Z=list(Zs),
#                    paraPen=list(Z1=list(S1=S1)),
#                    paraCon = list(Z1=Con.p),
#                    intrcpt.bam = F,
#                    model="linear",postselection = FALSE,
#                    est_beta_method = "multiridge")
#      })
#      temp <- data.frame("Method"="ecpc_scam", "Time"=mem[[2]], "Peak_memory"=mem[[4]],
#                         "n"=dim(AllData_n[[i]]$Xctd)[1], p=dim(AllData_n[[i]]$Xctd)[2])
#      df_n <- rbind(df_n, temp)
#  
#      #glmnet ridge
#      mem <- peakRAM({
#        fit.glmnet <- glmnet::cv.glmnet(y=AllData_n[[i]]$Y,x=AllData_n[[i]]$Xctd,
#                                        family='gaussian',alpha=0)
#        beta.glmnet <- coef(fit.glmnet,s="lambda.min", exact=TRUE)
#      })
#      temp <- data.frame("Method"="glmnet_ridge", "Time"=mem[[2]], "Peak_memory"=mem[[4]],
#                         "n"=dim(AllData_n[[i]]$Xctd)[1], p=dim(AllData_n[[i]]$Xctd)[2])
#      df_n <- rbind(df_n, temp)
#  
#  
#      #glmnet lasso
#      mem <- peakRAM({
#        fit.glmnet <- glmnet::cv.glmnet(y=AllData_n[[i]]$Y,x=AllData_n[[i]]$Xctd,
#                                        family='gaussian',alpha=1)
#        beta.glmnet <- coef(fit.glmnet,s="lambda.min", exact=TRUE)
#      })
#      temp <- data.frame("Method"="glmnet_lasso", "Time"=mem[[2]], "Peak_memory"=mem[[4]],
#                         "n"=dim(AllData_n[[i]]$Xctd)[1], p=dim(AllData_n[[i]]$Xctd)[2])
#      df_n <- rbind(df_n, temp)
#  
#      save(df_n, file=fname)
#    }
#  
#  }

## ---- eval=FALSE, class.source = 'fold-hide'----------------------------------
#  #Run for n=200, different p and check time and peak memory----
#  if(run_p){
#    fname <- paste(pathResults, 'res_p.Rdata', sep='')
#    df_p = data.frame("Method"=c(), "Time"=c(), "Peak_memory"=c(),
#                      "n" = c(), "p"=c())
#  
#    for(i in 1:length(p_all)){
#      p = p_all[i]
#  
#      #Compute co-data for p
#      Z <- abs(AllData_p[[i]]$beta)
#      ZI <- cbind(rep(1,p),abs(AllData_p[[1]]$beta))
#      Zs <- createZforSplines(values=Z, G=20, bdeg=3)
#      S1 <- createS(orderPen=2, G=20) #create 2nd order difference penalty matrix (same for all co-data)
#      Con.p <- createCon(G=20, shape="positive")
#  
#      #ecpc linear co-data model
#      mem <- peakRAM({
#        fit <- ecpc(Y=AllData_p[[i]]$Y,X=AllData_p[[i]]$Xctd,
#                    Z=list(ZI),
#                    model="linear",postselection=FALSE,est_beta_method = "multiridge")
#      })
#      temp <- data.frame("Method"="ecpc_linear", "Time"=mem[[2]], "Peak_memory"=mem[[4]],
#                         "n"=dim(AllData_p[[i]]$Xctd)[1], p=dim(AllData_p[[i]]$Xctd)[2])
#      df_p <- rbind(df_p, temp)
#  
#      #ecpc GAM co-data model
#      mem <- peakRAM({
#        fit <- ecpc(Y=AllData_p[[i]]$Y,X=AllData_p[[i]]$Xctd,
#                    Z=list(Zs),
#                    paraPen=list(Z1=list(S1=S1)),
#                    bam.method="ML", intrcpt.bam = F,
#                    model="linear",postselection=FALSE,
#                    est_beta_method = "multiridge")
#      })
#      temp <- data.frame("Method"="ecpc_gam", "Time"=mem[[2]], "Peak_memory"=mem[[4]],
#                         "n"=dim(AllData_p[[i]]$Xctd)[1], p=dim(AllData_p[[i]]$Xctd)[2])
#      df_p <- rbind(df_p, temp)
#  
#      #ecpc SCAM co-data model
#      mem <- peakRAM({
#        fit <- ecpc(Y=AllData_p[[i]]$Y,X=AllData_p[[i]]$Xctd,
#                    Z=list(Zs),
#                    paraPen=list(Z1=list(S1=S1)),
#                    paraCon = list(Z1=Con.p),
#                    intrcpt.bam = F,
#                    model="linear",postselection = FALSE,
#                    est_beta_method = "multiridge")
#      })
#      temp <- data.frame("Method"="ecpc_scam", "Time"=mem[[2]], "Peak_memory"=mem[[4]],
#                         "n"=dim(AllData_p[[i]]$Xctd)[1], p=dim(AllData_p[[i]]$Xctd)[2])
#      df_p <- rbind(df_p, temp)
#  
#      #glmnet ridge
#      mem <- peakRAM({
#        fit.glmnet <- glmnet::cv.glmnet(y=AllData_p[[i]]$Y,x=AllData_p[[i]]$Xctd,
#                                        family='gaussian',alpha=0)
#        beta.glmnet <- coef(fit.glmnet,s="lambda.min", exact=TRUE)
#      })
#      temp <- data.frame("Method"="glmnet_ridge", "Time"=mem[[2]], "Peak_memory"=mem[[4]],
#                         "n"=dim(AllData_p[[i]]$Xctd)[1], p=dim(AllData_p[[i]]$Xctd)[2])
#      df_p <- rbind(df_p, temp)
#  
#  
#      #glmnet lasso
#      mem <- peakRAM({
#        fit.glmnet <- glmnet::cv.glmnet(y=AllData_p[[i]]$Y,x=AllData_p[[i]]$Xctd,
#                                        family='gaussian',alpha=1)
#        beta.glmnet <- coef(fit.glmnet,s="lambda.min", exact=TRUE)
#      })
#      temp <- data.frame("Method"="glmnet_lasso", "Time"=mem[[2]], "Peak_memory"=mem[[4]],
#                         "n"=dim(AllData_p[[i]]$Xctd)[1], p=dim(AllData_p[[i]]$Xctd)[2])
#      df_p <- rbind(df_p, temp)
#  
#      save(df_p, file=fname)
#    }
#  
#  }

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#Load data for plots----
#different n, p=5000 fixed
fname <- paste(pathResults, 'res_n.Rdata', sep='')
load(fname)
df_n$Method <- factor(df_n$Method,levels=unique(df_n$Method),
                        labels=unique(df_n$Method))

#different p, n=200 fixed
fname <- paste(pathResults, 'res_p.Rdata', sep='')
load(fname)
df_p$Method <- factor(df_p$Method,levels=unique(df_p$Method),
                      labels=unique(df_p$Method))

## ---- fig.width=8, fig.height=5-----------------------------------------------
p1 <- ggplot(df_n)+
  aes(x=n,y=Time,col=Method,linetype=Method,shape=Method)+
  geom_line(linewidth=1.2)+
  geom_point(size=ps*1.5)+
  scale_color_manual(values=c('black','black','black','grey40','grey40'))+
  scale_linetype_manual(values=c(1,2,3,1,2))+
  scale_shape_manual(values=c(15,16,17,15,16))+
  #scale_y_log10()+
  labs(y="Time (s)",x="n")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts),
        legend.key.width = unit(2, "cm"))#,
p1

## ---- fig.width=8, fig.height=5-----------------------------------------------
p2 <- ggplot(df_n)+
  aes(x=n,y=Peak_memory,col=Method,linetype=Method,shape=Method)+
  geom_line(linewidth=1.2)+
  geom_point(size=ps*1.5)+
  scale_color_manual(values=c('black','black','black','grey40','grey40'))+
  scale_linetype_manual(values=c(1,2,3,1,2))+
  scale_shape_manual(values=c(15,16,17,15,16))+
  #scale_y_log10()+
  labs(y="Peak memory (MiB)",x="n")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts),
        legend.key.width = unit(2, "cm"))#,
p2


## ---- fig.width=8, fig.height=5-----------------------------------------------
p3 <- ggplot(df_p)+
  aes(x=p,y=Time,col=Method,linetype=Method,shape=Method)+
  geom_line(linewidth=1.2)+
  geom_point(size=ps*1.5)+
  scale_color_manual(values=c('black','black','black','grey40','grey40'))+
  scale_linetype_manual(values=c(1,2,3,1,2))+
  scale_shape_manual(values=c(15,16,17,15,16))+
  #scale_y_log10()+
  labs(y="Time (s)",x="p")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts),
        legend.key.width = unit(2, "cm"))#,
p3

## ---- fig.width=8, fig.height=5-----------------------------------------------
p4 <- ggplot(df_p)+
  aes(x=p,y=Peak_memory,col=Method,linetype=Method,shape=Method)+
  geom_line(linewidth=1.2)+
  geom_point(size=ps*1.5)+
  scale_color_manual(values=c('black','black','black','grey40','grey40'))+
  scale_linetype_manual(values=c(1,2,3,1,2))+
  scale_shape_manual(values=c(15,16,17,15,16))+
  #scale_y_log10()+
  labs(y="Peak memory (MiB)",x="p")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts),
        legend.key.width = unit(2, "cm"))#,
p4

## ---- fig.width=12, fig.height=8----------------------------------------------
p <- ggarrange(
  p1, p2, p3, p4,
  common.legend = TRUE, legend = "bottom"
)
p

