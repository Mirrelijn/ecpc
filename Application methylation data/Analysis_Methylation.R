##Verlaat methylation data (n=64, p=2720)
#Two Groupings: 
#  1. CpG islands (ridge hypershrinkage)
#  2. P-values from external study (hierarchical lasso + ridge hypershrinkage)
#NOTE 1: in Rstudio; press ALT+o for overview of code section headings
#NOTE 2: if whole file is run, analysis (on CV folds&subsamples for several methods) is skipped,
#      but results are loaded (if saved previously) and plots are saved to the working directory

#Load libraries----
#setwd("...") #directory in which results and figures are saved
pathData <- "" #path to data and co-data
library(ecpc)
pathResults <- "./Results/" #if analyses have been run previously, where can we find the results?
# load libraries required for ecpc
library(MASS)
library(penalized)
library(glmnet)
library(mvtnorm)
library(gglasso)
library(mgcv)
library(CVXR)
library(GRridge)
library(randomForest)
library(expm)
library(Rsolnp)
library(foreach)
library(doParallel)
library(graper)

#load libraries needed for storing and plotting results
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggraph)
library(igraph)
library(tidyverse)

#optional: 
if(0){ #set to 1 to setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(5) #not to overload your computer
  registerDoParallel(cl)
}

#Load data----
load(paste(pathData,"inputWridgeNonNAFU_UPDOWN.Rdata",sep=""))

dim(betasfcen) #centered methylation data
sds <- apply(betasfcen,1,sd) #compute standard deviations
betasfstand <- betasfcen/c(sds) #pxn: 2720 x 64. Standardize data
pvalFarkas <- selsumm[,5] #p-values from external study needed for co-data grouping

n <- dim(betasfcen)[2] #number of samples
p <- dim(betasfcen)[1] #number of covariates

Y <- fl #response: numeric values of 0/1
Xstd<-t(betasfstand) #standardised methylation data

#Define co-data groupings and types of hypershrinkage----
#ecpc takes as input the variable groupings: a list with m nested lists for m co-data groupings;
# each of those m nested lists contain the G_m groups of one grouping;
# each of those G_m groups contain the indices of covariates that belong to that group

#Grouping 1: based on distance to CpG-islands
load(paste(pathData,"Annotation.Rdata",sep="")) #load annotation data
CpGann <- factor(as.character(Annotation[whfark0,10])) #create a factor of the CpG-island groups
GroupingCpG <- CreatePartition(CpGann) #use the CreatePartition function from GRridge to create a list with groups corresponding to the factor levels

#Grouping 2: based on external p-values
#ecpc requires a discretised version of the continuous co-data;
# we use an adaptive discretisation, by using hierarchical groups of varying grid size,
# and using hierarchical lasso shrinkage on group level to select hierarchical groups.
#First create a list with the groups of covariates varying in size;
# splitMedian splits continuous co-data recursively at the median to form two new groups, 
# split="left" splits only the lower half group
GroupingPvals<-splitMedian(values=pvalFarkas,index=1:p,minGroupSize=50,split="left")
#Then define the hierarchy by forming groups on group level
HierarchyPvals <- obtainHierarchy(GroupingPvals) 

#optional: define grouping for standard deviations that may be added
GroupingSds <- CreatePartition(sds,ngroup=10,uniform=T) 

#Define groupings used in dense setting (and combined with posterior selection; in sparse setting)
GroupingsAll <- list(CpGislands=GroupingCpG,
                     Pvalues = GroupingPvals)
groupings.grouplvl <- list(NULL,HierarchyPvals) #Define groups on group level, in this case used to define the hierarchy on the groups of the second grouping
hypershrinkage <- c("ridge","hierLasso,ridge") #define type of hypershrinkage: ridge for the first grouping, a combination of hierarchical lasso and ridge for the second

#group sparse setting:
# combine lasso on group level for selecting groups with ridge to estimate weights of selected groups
if(0){ #set to 1 if group sparse results are preferred
  GroupingsAll <- list(CpGislands=GroupingCpG,
                       Pvalues = GroupingPvals)
  groupings.grouplvl <- list(NULL,NULL) 
  hypershrinkage <- c("lasso,ridge","lasso,ridge")
  
  #or dense setting with sds added as third grouping
  GroupingsAll <- list(CpGislands=GroupingCpG,
                       Pvalues = GroupingPvals,
                       Sds= GroupingSds)
  groupings.grouplvl <- list(NULL,HierarchyPvals,NULL) 
  hypershrinkage <- c("ridge","hierLasso,ridge","ridge")
}

#Groupings for GRridge (called "partitions" in GRridge function): CpG-islands and leaf groups of the hierarchical p-values groups
#Grouping: only use leaves of the hierarchical tree of grouping 2
leaves <- which(sapply(obtainHierarchy(HierarchyPvals),length)==1) #find leaves (groups that have no children)
GroupingPvalsGRridge<-GroupingPvals[leaves] #create grouping with leaf groups only
ord<-order(sapply(GroupingPvalsGRridge,function(x){mean(pvalFarkas[x])})) #find order of the p-value groups
GroupingPvalsGRridge <- GroupingPvalsGRridge[ord] #sort groups in decreasing p-values
partitionsGR <- list("CpG"=GroupingCpG,"Pval"=GroupingPvalsGRridge) 
monotoneGR <- c(FALSE,TRUE) #in GRridge, monotonicity can be imposed on the grouping weights

#annotation for graper 
#(similar to grouping used in ecpc, but slightly different format, for non-overlapping groups only)
#p-values grouping has been shown using ecpc to be informative -> use leaves of that grouping
if(1){
  #use leaves of p-values only as graper cannot handle overlapping groups nor multiple groupings
  annotGraper <- rep(0,p)
  for(i in 1:length(GroupingPvalsGRridge)){
    annotGraper[GroupingPvalsGRridge[[i]]] <- i 
  }
}else{
  values <- pvalFarkas
  GroupingPvals.10groups <- CreatePartition(values,ngroup=10,uniform=T) #10 groups
  annotGraper <- rep(0,p)
  for(i in 1:length(GroupingPvals.10groups)){
    annotGraper[GroupingPvals.10groups[[i]]] <- i 
  }
}


#Define folds and variables for the CV----
logname<-paste("logCVVerlaat.txt",sep="")
write(paste(Sys.time(),logname),file=logname,append=T)

nfolds<-20
set.seed(43875492) #set.seed(8430)
folds2<-produceFolds(n,nfolds,Y,balance=T) #produce folds balanced in response
#save(folds2,file="foldsVerlaat3")
load("foldsVerlaat")
#lambdas <- sapply(Res,function(x) x$lambdaridge)
load("lambdasVerlaat") #optional: use lambdas saved after run ecpc for fair comparison ecpc with different hypershrinkage
maxSel <- c(2:29,seq(30,100,by=10)) #define a range of maximum number of selected variables for posterior selection

#Define subsamples and variables for subsamples----
maxSel <- 25 #select 25 or 50 variables
alp <- c(0.3,0.8) #alpha variable used in elastic net 

set.seed(47198)
nSmpl <- 50 #number of subsample iterations
ind0 <- which(Y==0)
ind1 <- which(Y==1)
sizeSample0 <- floor(length(ind0)*2/3)
sizeSample1 <- floor(length(ind1)*2/3)
subsamples0 <- replicate(nSmpl,sample(ind0,size=sizeSample0,replace=F))
subsamples1 <- replicate(nSmpl,sample(ind1,size=sizeSample1,replace=F))
#subsamples <- rbind(subsamples0,subsamples1) #stratified subsamples
#save(subsamples,file="subsamples")
load("subsamplesVerlaat")

#Do CV for ecpc and ordinary ridge----
fname <- "CVMethylation_ecpcA" #dense setting and sparse setting
#fname <- "CVMethylation_ecpcB" #group sparse setting

if(0){ #set to 1 to perform analysis or run code inside block manually
  grpsno <- c(unlist(sapply(GroupingsAll,function(x){1:length(x)}))) #vector with group numbers in all groupings
  grpngsno <- c(unlist(sapply(1:length(GroupingsAll),function(i){rep(i,length(GroupingsAll[[i]]))}))) #vector with groupings numbers
  dfGrps<-data.frame() #data frame in which group and grouping weights are stored
  dfBeta <- data.frame()
  df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
  Res<-list() #list in which raw output of ecpc is stored (e.g. estimated regression coefficients)
  #for(i in 1:nfolds){
  finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                         .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                       "CVXR","GRridge","expm","Rsolnp","ecpc")) %dopar% {
    tic<-proc.time()[[3]]
    Res[[i]]<-ecpc(Y[-folds2[[i]]],Xstd[-folds2[[i]],],GroupingsAll,hypershrinkage=hypershrinkage,
                   Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],#lambda=lambdas[i],
                   groupings.grouplvl = groupings.grouplvl,postselection="elnet+dense",maxsel=maxSel,normalise=T)
    Res[[i]]$timeGR <- proc.time()[[3]]-tic
    df2<-data.frame("Ypred"=c(Res[[i]]$Ypred,Res[[i]]$Ypredridge,c(Res[[i]]$YpredPost)))
    df2$Method <- rep(c("ecpc","ordinary.ridge",paste("ecpc",maxSel,"vars",sep="")),each=length(folds2[[i]]))
    df2$NumberSelectedVars <- rep(c(sum(Res[[i]]$beta!=0),p,maxSel),each=length(folds2[[i]]))
    df2$Fold <- i
    df2$Sample <- rep(folds2[[i]],2+length(maxSel))
    df2$Time <-  Res[[i]]$timeGR
    df2$Truth <- rep(Y[folds2[[i]]],2+length(maxSel))
    #df<-rbind(df,df2)

    df3<-data.frame("Group"=grpsno,
                    "Grouping"=grpngsno,
                    "Group weight"=Res[[i]]$gamma,
                    "Grouping weight"=Res[[i]]$w[grpngsno])
    df3$Tau.ecpc <- Res[[i]]$tauglobal #global tau^2
    df3$Tau.ridge <- 1/Res[[i]]$lambdaridge #ordinary ridge tau^2
    df3$Method <- "ecpc"
    df3$Fold <- i
    #dfGrps<-rbind(dfGrps,df3)
    
    whichPost2 <- apply(Res[[i]]$betaPost,2,function(x){which(x!=0)})
    betaPost2 <- unlist(lapply(1:dim(Res[[i]]$betaPost)[2],function(x){Res[[i]]$betaPost[whichPost2[[x]],x]}))
    setting2 <- lapply(1:length(whichPost2),function(x){rep(x,length(whichPost2[[x]]))})
    MaxSel2 <- lapply(1:length(whichPost2),function(x){rep(length(whichPost2[[x]]),length(whichPost2[[x]]))})
    dfBeta2<-data.frame("Betas"=c(betaPost2),
                        "whichPost"=c(unlist(whichPost2)),
                        "Setting"=c(unlist(setting2)),
                        "MaxSel"=c(unlist(MaxSel2)),
                        "Method"=c(rep("ecpc+selection",length(betaPost2))))

    dfBeta2$Fold <- i
    #dfBeta<-rbind(dfBeta,dfBeta2)
    
    write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
    
    list("Res"=Res,"df"=df2,"dfGrps"=df3,"dfBeta"=dfBeta2)
  }
  
  Res <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]][[i]])
  df2 <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]])
  dfGrps2 <- lapply(1:nfolds,function(i) finalMatrix[i,3][[1]])
  dfBeta2 <- lapply(1:nfolds,function(i) finalMatrix[i,4][[1]])
  df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
  dfGrps <- dfGrps2[[1]]; for(i in 2:nfolds) dfGrps <- rbind(dfGrps,dfGrps2[[i]])
  dfBeta <- dfBeta2[[1]]; for(i in 2:nfolds) dfBeta <- rbind(dfBeta,dfBeta2[[i]])
  save(Res,df,dfGrps,dfBeta,file="CVVerlaat2")
  stopCluster(cl); rm(cl)
  
  #data frame with summary statistics
  Summdf <- df %>% group_by(Method) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                  CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                  NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  df$Method<-as.factor(df$Method)
  dfGrps$Group<-as.factor(dfGrps$Group)
  
  dfROC<-data.frame()
  for(i in levels(df$Method)){
    temp<-data.frame()
    cutoffs<-rev(seq(0,1,by=0.001))
    rocGR <- roc(probs=df$Ypred[df$Method==i],true=df$Truth[df$Method==i],cutoffs=cutoffs)
    temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
    temp$Method <- i
    temp$AUC<-c(auc(rocGR))
    temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
    dfROC<-rbind(dfROC,temp)
  }
  dfAUC <- dfROC %>% group_by(Method) %>% summarise(AUC=mean(AUC),
                                                    NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  #plot(dfAUC$NumberSelectedVars,dfAUC$AUC)
  
  if(all(Summdf$NumberSelectedVars==dfAUC$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
    Summdf$AUC <- dfAUC$AUC
  }
  dfGrps$Group<-as.factor(dfGrps$Group)
  save(df,dfGrps,Summdf,dfROC,dfAUC,Res,dfBeta,file=fname) #this comparison 
}


#Do CV for elastic net-----
alp<-c(0.3,0.8)
fname <- "CVMethylation_ElNet" #file in which the results are stored

if(0){
  dfElNet <- data.frame()
  dfBetaElNet <- data.frame()
  ResElNet <- list()
  for(j in 1:length(alp)){
    for(i in 1:nfolds){
      ResElNet[[i]]<-glmnet(Xstd[-folds2[[i]],],Y[-folds2[[i]]],alpha=alp[j],family="binomial",
                       standardize = F,intercept=T) 
      Ypred<-predict(ResElNet[[i]],newx=Xstd[folds2[[i]],], type="response")
      
      df2<-data.frame("Ypred"=c(Ypred))
      df2$Method <- "ElNet"
      df2$NumberSelectedVars <- rep(c(ResElNet[[i]]$df),each=length(folds2[[i]]))
      df2$Setting <- rep(1:100,each=length(folds2[[i]]))
      df2$Fold <- i
      df2$Sample <- rep(folds2[[i]],100)
      df2$Time <-  Res[[i]]$timeGR
      df2$Truth <- rep(Y[folds2[[i]]],100)
      df2$Alpha <- alp[j]
      dfElNet<-rbind(dfElNet,df2)
      
      whichPost2<-apply(ResElNet[[i]]$beta,2,function(x){which(x!=0)})
      betaPost2 <- unlist(lapply(1:length(whichPost2),function(x){ResElNet[[i]]$beta[whichPost2[[x]],x]}))
      setting2 <- lapply(1:length(whichPost2),function(x){rep(x,length(whichPost2[[x]]))})
      MaxSel2 <- lapply(1:length(whichPost2),function(x){rep(length(whichPost2[[x]]),length(whichPost2[[x]]))})
      dfBetaElNet2<-data.frame("Betas"=c(betaPost2),
                               "whichPost"=c(unlist(whichPost2)),
                               "Setting"=c(unlist(setting2)),
                               "MaxSel"=c(unlist(MaxSel2)),
                               "Method"=c(rep("ElNet",length(betaPost2))),
                               "ShrinkMethod"=c(rep("ElNet",length(betaPost2))),
                               "SelectMethod"=c(rep("ElNet",length(betaPost2)))
      )
      dfBetaElNet2$Fold <- i
      dfBetaElNet2$Alpha <- alp[j]
      dfBetaElNet<-rbind(dfBetaElNet,dfBetaElNet2)
      
      write(paste(Sys.time(),"fold",i,"of",nfolds,"done, elastic net"),file=logname,append=T)
      
    }
  }
  
  #data frame with summary statistics
  SummdfElNet <- dfElNet %>% group_by(Alpha,Setting,Method) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                                          CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                                          NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  dfElNet$Method<-as.factor(dfElNet$Method)
  
  dfROCElNet<-data.frame()
  for(i in levels(dfElNet$Method)){
    for(j in unique(dfElNet$Alpha)){
      for(k in unique(dfElNet$Setting)){
        temp<-data.frame()
        cutoffs<-rev(seq(0,1,by=0.001))
        rocGR <- roc(probs=dfElNet$Ypred[dfElNet$Method==i & dfElNet$Alpha==j & dfElNet$Setting==k],
                     true=dfElNet$Truth[dfElNet$Method==i & dfElNet$Alpha==j & dfElNet$Setting==k],
                     cutoffs=cutoffs)
        temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
        temp$Method <- i
        temp$Alpha <- j
        temp$Setting <- k
        temp$AUC<-c(auc(rocGR))
        temp$NumberSelectedVars<-mean(dfElNet$NumberSelectedVars[dfElNet$Method==i & dfElNet$Alpha==j & dfElNet$Setting==k])
        dfROCElNet<-rbind(dfROCElNet,temp)
      }
    }
  }
  dfAUCElNet <- dfROCElNet %>% group_by(Method,Alpha,Setting) %>% 
    summarise(AUC=mean(AUC),NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  #plot(dfAUCElNet$NumberSelectedVars,dfAUCElNet$AUC)
  
  if(all(SummdfElNet$NumberSelectedVars==dfAUCElNet$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
    SummdfElNet$AUC <- dfAUCElNet$AUC
  }
  
  save(dfElNet,SummdfElNet,dfROCElNet,dfAUCElNet,ResElNet,dfBetaElNet,file=fname) 
}


#Do CV for GRridge----
fname <- "CVMethylation_GRridge" #file in which the results are stored

if(0){
  grpsnoGR <- unlist(sapply(partitionsGR,function(x){1:length(x)}))
  grpngsnoGR <- unlist(sapply(1:length(partitionsGR),function(i){rep(i,length(partitionsGR[[i]]))}))
  dfGrpsGRridge<-data.frame()
  dfGRridge<-data.frame()
  dfBetaGRridge <- data.frame()
  ResGRridge<-list()
  for(i in 1:nfolds){
    tic<-proc.time()[[3]]
    ResGRridge[[i]] <- grridge(t(Xstd[-folds2[[i]],]),Y[-folds2[[i]]], partitions = partitionsGR,
                               monotone=monotoneGR ,optl=lambdas[i])
    time <- proc.time()[[3]]-tic
    Xb<-Xstd[folds2[[i]],]%*%ResGRridge[[i]]$betas+ResGRridge[[i]]$predobj$GroupRegul@unpenalized
    Ypredgrridge<-1/(1+exp(-Xb)) #predictions
    
    df2<-data.frame("Ypred"=Ypredgrridge)
    df2$Method <- "GRridge"
    df2$NumberSelectedVars <- p
    df2$Fold <- i
    df2$Sample <- folds2[[i]]
    df2$Time <-  time
    df2$Truth <- Y[folds2[[i]]]
    df2$Setting <- p
    dfGRridge<-rbind(dfGRridge,df2)
    
    df3<-data.frame("Group"=grpsnoGR,
                    "Grouping"=grpngsnoGR,
                    "Group weight"=1/unlist(ResGRridge[[i]]$lambdamults),
                    "Grouping weight"=rep(1,length(grpngsnoGR)))
    df3$Tau.GR <- 1/ResGRridge[[i]]$optl #global tau
    df3$Method <- "GRridge"
    df3$Fold <- i
    dfGrpsGRridge<-rbind(dfGrpsGRridge,df3)
    
    #GRridge post-selection 
    penalties <- ResGRridge[[i]]$lambdamultvec[,dim(ResGRridge[[i]]$lambdamultvec)[2]]*ResGRridge[[i]]$optl
    post2 <- postSelect(X=Xstd[-folds2[[i]],],Y=Y[-folds2[[i]]],
                           beta=ResGRridge[[i]]$betas,intrcpt=ResGRridge[[i]]$predobj$GroupRegul@unpenalized,
                           postselection="elnet+dense",maxsel=maxSel, 
                           penalties=penalties,model="logistic",
                           tauglobal=1/ResGRridge[[i]]$optl,
                           X2=Xstd[folds2[[i]],],Y2=Y[folds2[[i]]])
    df2<-data.frame("Ypred"=c(c(post2$YpredPost)))
    df2$Method<-rep(paste("GRridge",maxSel,"vars",sep=""),each=length(folds2[[i]]))
    df2$NumberSelectedVars <- c(rep(apply(post2$betaPost,2,function(x){sum(x!=0)}),each=length(folds2[[i]])))
    df2$Fold <- i
    df2$Sample <- folds2[[i]]
    df2$Time <-  time
    df2$Truth <- Y[folds2[[i]]]
    df2$Setting <- as.factor(rep(1:length(maxSel),each=length(folds2[[i]])))
    dfGRridge<-rbind(dfGRridge,df2)
    
    whichPost2 <- apply(post2$betaPost,2,function(x){which(x!=0)})
    betaPost2 <- unlist(lapply(1:dim(post2$betaPost)[2],function(x){post2$betaPost[whichPost2[[x]],x]}))
    setting2 <- lapply(1:length(whichPost2),function(x){rep(x,length(whichPost2[[x]]))})
    MaxSel2 <- lapply(1:length(whichPost2),function(x){rep(length(whichPost2[[x]]),length(whichPost2[[x]]))})
    dfBeta2<-data.frame("Betas"=c(betaPost2),
                        "whichPost"=c(unlist(whichPost2)),
                        "Setting"=c(unlist(setting2)),
                        "MaxSel"=c(unlist(MaxSel2)),
                        "Method"=c(rep("GRridge+postElNet",length(betaPost2))))
    dfBeta2$Fold <- i
    dfBetaGRridge<-rbind(dfBetaGRridge,dfBeta2)
    
    
    write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
  }
  
  SummdfGRridge <- dfGRridge %>% group_by(Method) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                  CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                  NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  dfGRridge$Method<-as.factor(dfGRridge$Method)
  dfGrpsGRridge$Group<- as.factor(dfGrpsGRridge$Group)
  
  dfROCGRridge<-data.frame()
  for(i in levels(dfGRridge$Method)){
    temp<-data.frame()
    cutoffs<-rev(seq(0,1,by=0.001))
    rocGR <- roc(probs=dfGRridge$Ypred[dfGRridge$Method==i],
                 true=dfGRridge$Truth[dfGRridge$Method==i],
                 cutoffs=cutoffs)
    temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
    temp$Method <- i
    temp$AUC<-c(auc(rocGR))
    temp$NumberSelectedVars<-mean(dfGRridge$NumberSelectedVars[dfGRridge$Method==i])
    dfROCGRridge<-rbind(dfROCGRridge,temp)
  }
  dfAUCGRridge <- dfROCGRridge %>% group_by(Method) %>% summarise(AUC=mean(AUC),
                                                    NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  #plot(dfAUCGRridge$NumberSelectedVars,dfAUCGRridge$AUC)
  
  if(all(SummdfGRridge$NumberSelectedVars==dfAUCGRridge$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
    SummdfGRridge$AUC <- dfAUCGRridge$AUC
  }
  save(dfGRridge,dfGrpsGRridge,SummdfGRridge,dfROCGRridge,dfAUCGRridge,ResGRridge,dfBetaGRridge,file=fname) #this comparison 
}



#Do CV for hierarchical lasso on covariate level----
fname <- "CVMethylation_HierLasso" #CpG islands and all hierarchical p-value groups
groupingHL <- c(GroupingCpG,GroupingPvals) 

if(0){
  grpsnoHL <- unlist(sapply(groupingsHL,function(x){1:length(x)}))
  grpngsnoHL <- unlist(sapply(1:length(groupingsHL),function(i){rep(i,length(groupingsHL[[i]]))}))
  dfGrpsHL<-data.frame()
  dfHL<-data.frame()
  ResHL<-list()
  for(i in 1:nfolds){
    tic<-proc.time()[[3]]
    ResHL[[i]] <- hierarchicalLasso(X=Xstd[-folds2[[i]],], Y=Y[-folds2[[i]]] , grouping=GroupingHL)
    ResHL[[i]]$time <- proc.time()[[3]]-tic
    Xb <- Xstd[folds2[[i]],]%*%ResHL[[i]]$betas + ResHL[[i]]$a0
    Ypred <- 1/(1+exp(-Xb))
    
    df2<-data.frame("Ypred"=Ypred)
    df2$Method <- "HierLasso"
    df2$NumberSelectedVars <- rep(sum(ResHL[[i]]$betas!=0),each=length(folds2[[i]]))
    df2$Fold <- i
    df2$Sample <- folds2[[i]]
    df2$Time <-  ResHL[[i]]$time
    df2$Truth <- Y[folds2[[i]]]
    dfHL<-rbind(dfHL,df2)
    
    df3<-data.frame("Group"=c(grpsnoHL),
                    "Grouping"=grpngsnoHL,
                    "Group weight"=ResHL[[i]]$group.weights,
                    "Grouping weight"=rep(1,length(grpngsnoHL)))
    df3$Tau.HL <- 1/ResHL[[i]]$lambda #overall tau
    df3$Method <- "HierLasso"
    df3$Fold <- i
    dfGrpsHL<-rbind(dfGrpsHL,df3)
    
    write(paste(Sys.time(),"fold",i,"of",nfolds,"done, Hierarchical Lasso"),file=logname,append=T)
    
  }
  
  SummdfHL <- dfHL %>% group_by(Method) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                       CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                       NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup
  
  dfHL$Method<-as.factor(dfHL$Method)
  dfGrpsHL$Group<-as.factor(dfGrpsHL$Group)
  
  dfROCHL<-data.frame()
  for(i in levels(dfHL$Method)){
    temp<-data.frame()
    cutoffs<-rev(seq(0,1,by=0.001))
    rocGR <- roc(probs=dfHL$Ypred[dfHL$Method==i],true=dfHL$Truth[dfHL$Method==i],cutoffs=cutoffs)
    temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
    temp$Method <- i
    temp$AUC<-c(auc(rocGR))
    temp$NumberSelectedVars<-mean(dfHL$NumberSelectedVars[dfHL$Method==i])
    dfROCHL<-rbind(dfROCHL,temp)
  }
  #dfAUCHL <- dfROCHL %>% group_by(Method) %>% summarise(AUC=mean(AUC),NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  dfAUCHL <- dfROCHL %>% group_by(Method) %>% summarise(AUC=mean(AUC),
                                                         NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  dfAUCHL<-dfHL %>% group_by(Method) %>% distinct(Fold, .keep_all=T) %>% select(NumberSelectedVars,Method) %>% 
    mutate(AUC=dfAUCHL$AUC[dfAUCHL$Method==unique(Method)]) %>% ungroup()
  save(ResHL,dfHL,dfGrpsHL,SummdfHL,dfROCHL,dfAUCHL,file=fname) #this comparison
}



#Do CV for group lasso on covariate level----
#use LOG penalty for overlapping groups
fname <- "CVMethylation_GL" #both p-value and cpg groups
GroupingGL <- c(GroupingCpG,GroupingPvalsGRridge)

if(0){
  grpsnoGL <- c(unlist(sapply(list(GroupingGL),function(x){1:length(x)})))
  grpngsnoGL <- c(unlist(sapply(1:length(list(GroupingGL)),function(i){rep(i,length(list(GroupingGL)[[i]]))})))
  dfGrpsGL<-data.frame()
  dfGL<-data.frame()
  ResGL<-list()
  for(i in 1:nfolds){
    tic<-proc.time()[[3]]
    ResGL[[i]] <- hierarchicalLasso(X=Xstd[-folds2[[i]],], Y=Y[-folds2[[i]]] , grouping=GroupingGL)
    ResGL[[i]]$time <- proc.time()[[3]]-tic
    Xb <- Xstd[folds2[[i]],]%*%ResGL[[i]]$betas + ResGL[[i]]$a0
    Ypred <- 1/(1+exp(-Xb))
    
    df2<-data.frame("Ypred"=Ypred)
    df2$Method <- "GroupLasso"
    df2$NumberSelectedVars <- rep(sum(ResGL[[i]]$betas!=0),each=length(folds2[[i]]))
    df2$Fold <- i
    df2$Sample <- folds2[[i]]
    df2$Time <-  ResGL[[i]]$time
    df2$Truth <- Y[folds2[[i]]]
    dfGL<-rbind(dfGL,df2)
    
    df3<-data.frame("Group"=c(grpsnoGL),
                    "Grouping"=grpngsnoGL,
                    "Group weight"=ResGL[[i]]$group.weights,
                    "Grouping weight"=rep(1,length(grpngsnoGL)))
    df3$Tau.GL <- 1/ResGL[[i]]$lambda #overall tau
    df3$Method <- "GroupLasso"
    df3$Fold <- i
    dfGrpsGL<-rbind(dfGrpsGL,df3)
    
    write(paste(Sys.time(),"fold",i,"of",nfolds,"done, Group Lasso"),file=logname,append=T)
    
  }
  
  SummdfGL <- dfGL %>% group_by(Method) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                      CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                      NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup
  dfGL$Method<-as.factor(dfGL$Method)
  dfGrpsGL$Group<-as.factor(dfGrpsGL$Group)
  
  dfROCGL<-data.frame()
  for(i in levels(dfGL$Method)){
    temp<-data.frame()
    cutoffs<-rev(seq(0,1,by=0.001))
    rocGR <- roc(probs=dfGL$Ypred[dfGL$Method==i],true=dfGL$Truth[dfGL$Method==i],cutoffs=cutoffs)
    temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
    temp$Method <- i
    temp$AUC<-c(auc(rocGR))
    temp$NumberSelectedVars<-mean(dfGL$NumberSelectedVars[dfGL$Method==i])
    dfROCGL<-rbind(dfROCGL,temp)
  }
  #dfAUCGL <- dfROCGL %>% group_by(Method) %>% summarise(AUC=mean(AUC),NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  dfAUCGL <- dfROCGL %>% group_by(Method) %>% summarise(AUC=mean(AUC),
                                                        NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  dfAUCGL<-dfGL %>% group_by(Method) %>% distinct(Fold, .keep_all=T) %>% select(NumberSelectedVars,Method) %>%
    mutate(AUC=dfAUCGL$AUC[dfAUCGL$Method==unique(Method)]) %>% ungroup()
  save(ResGL,dfGL,dfGrpsGL,SummdfGL,dfROCGL,dfAUCGL,file=fname) #this comparison
}

#Do CV for random forest----
fnameRF <- "CVMethylation_RF" #file in which the results are stored"

if(0){
  dfRF<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
  ResRF<-list() #list in which raw output of ecpc is stored (e.g. estimated regression coefficients)
  YF <- factor(Y)
  for(i in 1:nfolds){
    tic<-proc.time()[[3]]
    ResRF[[i]]<-randomForest(y=YF[-folds2[[i]]],x=Xstd[-folds2[[i]],],
                             ytest=YF[folds2[[i]]],xtest=Xstd[folds2[[i]],],keep.forest = T)
    Ypred <- apply(predict(ResRF[[i]],newdata=Xstd[folds2[[i]],],type="prob"),1,max)
    ResRF[[i]]$timeGR <- proc.time()[[3]]-tic
    df2<-data.frame("Ypred"=Ypred) 
    df2$Method <- "RandomForest"
    df2$NumberSelectedVars <- p
    df2$Fold <- i
    df2$Sample <- folds2[[i]]
    df2$Time <-  ResRF[[i]]$timeGR
    df2$Truth <- Y[folds2[[i]]]
    dfRF<-rbind(dfRF,df2)
    
    write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
  }
  
  #data frame with summary statistics
  SummdfRF <- dfRF %>% group_by(Method) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                      CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                      NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  dfRF$Method<-as.factor(dfRF$Method)
  
  dfROCRF<-data.frame()
  for(i in levels(dfRF$Method)){
    temp<-data.frame()
    cutoffs<-rev(seq(0,1,by=0.001))
    rocGR <- roc(probs=dfRF$Ypred[dfRF$Method==i],true=dfRF$Truth[dfRF$Method==i],cutoffs=cutoffs)
    temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
    temp$Method <- i
    temp$AUC<-c(auc(rocGR))
    temp$NumberSelectedVars<-mean(dfRF$NumberSelectedVars[dfRF$Method==i])
    dfROCRF<-rbind(dfROCRF,temp)
  }
  dfAUCRF <- dfROCRF %>% group_by(Method) %>% summarise(AUC=mean(AUC),NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  #plot(dfAUCRF$NumberSelectedVars,dfAUCRF$AUC)
  
  if(all(SummdfRF$NumberSelectedVars==dfAUCRF$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
    SummdfRF$AUC <- dfAUCRF$AUC
  }
  
  save(dfRF,SummdfRF,dfROCRF,dfAUCRF,ResRF,file=fnameRF) #this comparison 
}


#Do CV for graper----
fname <- "CVMethylation_graper" #leaf groups pvalues, sparse
#fname <- "CVMethylation_graperA" #leaf groups pvalues, dense (spikeslab=F)

if(0){
  grpsno <- 1:max(annotGraper)
  grpngsno <- rep(2,max(annotGraper))
  dfGrps<-data.frame()
  df<-data.frame()
  dfBeta <- data.frame()
  Res<-list()
  for(i in 1:nfolds){
    tic<-proc.time()[[3]]
    Res[[i]] <- graper(Xstd[-folds2[[i]],], Y[-folds2[[i]]], annotGraper, verbose=F,
                       family="binomial", calcELB=FALSE,standardize = F,spikeslab =F,factoriseQ = T)
    time <- proc.time()[[3]]-tic
    Xb<-Xstd[folds2[[i]],]%*%coef(Res[[i]], include_intercept=FALSE)+Res[[i]]$intercept
    Ypred<-1/(1+exp(-Xb)) #predictions
    #Ypred<- predict(Res[[i]], Xstd[folds2[[i]],]) #predictions
    
    df2<-data.frame("Ypred"=Ypred)
    df2$Method <- "graper"
    df2$NumberSelectedVars <- sum(coef(Res[[i]], include_intercept=FALSE)!=0)
    df2$Fold <- i
    df2$Sample <- folds2[[i]]
    df2$Time <-  time
    df2$Truth <- Y[folds2[[i]]]
    df2$Setting <- p
    df<-rbind(df,df2)
    
    df3<-data.frame("Group"=grpsno,
                    "Grouping"=grpngsno,
                    "Group weight"=1/unlist(Res[[i]]$EW_gamma),
                    "Grouping weight"=rep(1,length(grpngsno)))
    df3$Tau.graper <- 1 #overall tau
    df3$Method <- "graper"
    df3$Fold <- i
    dfGrps<-rbind(dfGrps,df3)
    
    
    #write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
  }
  
  
  Summdf <- df %>% group_by(Method) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                  CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                  NumberSelectedVars=mean(NumberSelectedVars))%>% ungroup()
  
  df$Method<-as.factor(df$Method)
  dfGrps$Group<- as.factor(dfGrps$Group)
  
  dfROC<-data.frame()
  for(i in levels(df$Method)){
    temp<-data.frame()
    cutoffs<-rev(seq(0,1,by=0.001))
    rocGR <- roc(probs=df$Ypred[df$Method==i],
                 true=df$Truth[df$Method==i],
                 cutoffs=cutoffs)
    temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
    temp$Method <- i
    temp$AUC<-c(auc(rocGR))
    temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
    dfROC<-rbind(dfROC,temp)
  }
  dfAUC <- dfROC %>% group_by(Method) %>% 
    summarise(AUC=mean(AUC),NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  #plot(dfAUC$NumberSelectedVars,dfAUC$AUC)
  
  if(all(Summdf$NumberSelectedVars==dfAUC$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
    Summdf$AUC <- dfAUC$AUC
  }
  save(df,dfGrps,Summdf,dfROC,dfAUC,Res,file=fname) #this comparison 
}



#Do posterior selection for ecpc on CV folds----
maxSel <- c(2:29,seq(30,100,by=10))
postselection <- "BRmarginal,dense"
postselection <- "BRmarginal,sparse"
postselection <- "elnet,dense"
postselection <- "elnet,sparse"
postselection <- "DSS"

if(0){
  fname <- paste("CVMethylation_ecpcA+",postselection,sep="")
  load(paste(pathResults,"CVMethylation_ecpcA",sep=""))
  dfPost <- data.frame()
  dfBetaPost <- data.frame()
  for(i in 1:nfolds){
    post2 <- postSelect(X=Xstd[-folds2[[i]],],Y=Y[-folds2[[i]]],
                           beta=Res[[i]]$beta,intrcpt=Res[[i]]$intercept,
                           postselection=postselection,maxsel=maxSel, 
                           penalties=Res[[i]]$penalties,model="logistic",
                           tauglobal=Res[[i]]$tauglobal,
                           X2=Xstd[folds2[[i]],],Y2=Y[folds2[[i]]])
    df2<-data.frame("Ypred"=c(c(post2$YpredPost)))
    df2$Method<-rep(paste("ecpc+",postselection,maxSel,"vars",sep=""),each=length(folds2[[i]]))
    df2$NumberSelectedVars <- c(rep(apply(post2$betaPost,2,function(x){sum(x!=0)}),each=length(folds2[[i]])))
    df2$Fold <- i
    df2$Sample <- rep(folds2[[i]],length(maxSel))
    df2$Time <-  NaN
    df2$Truth <- rep(Y[folds2[[i]]],length(maxSel))
    df2$Setting <- as.factor(rep(1:length(maxSel),each=length(folds2[[i]])))
    dfPost<-rbind(dfPost,df2)
    
    whichPost2 <- apply(post2$betaPost,2,function(x){which(x!=0)})
    betaPost2 <- unlist(lapply(1:dim(post2$betaPost)[2],function(x){post2$betaPost[whichPost2[[x]],x]}))
    setting2 <- lapply(1:length(whichPost2),function(x){rep(x,length(whichPost2[[x]]))})
    MaxSel2 <- lapply(1:length(whichPost2),function(x){rep(length(whichPost2[[x]]),length(whichPost2[[x]]))})
    dfBeta2<-data.frame("Betas"=c(betaPost2),
                        "whichPost"=c(unlist(whichPost2)),
                        "Setting"=c(unlist(setting2)),
                        "MaxSel"=c(unlist(MaxSel2)),
                        "Method"=c(rep(paste("ecpc+",postselection,sep=""),length(betaPost2))))
    dfBeta2$Fold <- i
    dfBetaPost<-rbind(dfBetaPost,dfBeta2)
    
    write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
  }
  
  #data frame with summary statistics
  SummdfPost <- dfPost %>% group_by(Method) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                  CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                  NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  dfPost$Method<-as.factor(dfPost$Method)
  
  dfROCPost<-data.frame()
  for(i in levels(dfPost$Method)){
    temp<-data.frame()
    cutoffs<-rev(seq(0,1,by=0.001))
    rocGR <- roc(probs=dfPost$Ypred[dfPost$Method==i],true=dfPost$Truth[dfPost$Method==i],cutoffs=cutoffs)
    temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
    temp$Method <- i
    temp$AUC<-c(auc(rocGR))
    temp$NumberSelectedVars<-mean(dfPost$NumberSelectedVars[dfPost$Method==i])
    dfROCPost<-rbind(dfROCPost,temp)
  }
  dfAUCPost <- dfROCPost %>% group_by(Method) %>% summarise(AUC=mean(AUC),NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  #plot(dfAUCPost$NumberSelectedVars,dfAUCPost$AUC)
  
  if(all(SummdfPost$NumberSelectedVars==dfAUCPost$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
    SummdfPost$AUC <- dfAUCPost$AUC
  }
  save(dfPost,SummdfPost,dfROCPost,dfAUCPost,dfBetaPost,file=fname) #this comparison 
}



#Do subsampling to compare overlap in selected variable sets----
logname<-paste("logSubsamplingVerlaat.txt",sep="")
write(paste(Sys.time(),logname),file=logname,append=T)

if(0){
  Res<-list()
  which25ecpc<-list("ecpc"=list())
  grpsno <- unlist(sapply(GroupingsAll,function(x){1:length(x)})) #vector with group numbers in all groupings
  grpngsno <- unlist(sapply(1:length(GroupingsAll),function(i){rep(i,length(GroupingsAll[[i]]))})) #vector with grouping numbers
  dfGrps<-data.frame() #data frame in which group and grouping weights are stored
  for(i in 1:nSmpl){
    #ecpc
    tic<-proc.time()[[3]]
    Res[[i]]<-ecpc(Y[subsamples[,i]],Xstd[subsamples[,i],],GroupingsAll,hypershrinkage=hypershrinkage,
                   groupings.grouplvl=groupings.grouplvl,postselection="elnet+dense",maxsel=maxSel)
    Res[[i]]$time<- proc.time()[[3]]-tic
    which25ecpc$ecpc[[i]]<-apply(Res[[i]]$betaPost,2,function(x){which(x!=0)})
    write(paste(Sys.time(),"Subsample ",i," done, ecpc",sep=""),file=logname,append=T)
    
    df3<-data.frame("Group"=grpsno,
                    "Grouping"=grpngsno,
                    "Group weight"=Res[[i]]$gamma,
                    "Grouping weight"=Res[[i]]$w[grpngsno])
    df3$Tau.ecpc <- Res[[i]]$tauglobal #global tau
    df3$Tau.ridge <- 1/Res[[i]]$lambdaridge #ordinary ridge tau
    df3$Method <- "ecpc"
    df3$Fold <- i
    dfGrps<-rbind(dfGrps,df3)
    
    save(which25ecpc,Res,dfGrps,file="which25ecpc")
  }
  
  which25EN <-list()
  ResElnet <- list()
  for(i in 1:length(alp)){
    which25EN[[i]]<-list()
    ResElnet[[i]] <- list()
    for(j in 1:nSmpl){
      ResElnet[[i]][[j]]<-list()
      
      tic<-proc.time()[[3]]
      ResElnet[[i]][[j]]$fit<-glmnet(Xstd[subsamples[,j],],Y[subsamples[,j]],alpha=alp[i],family="binomial",
                                     standardize = F,intercept=T) 
      
      #then find lambda values that correspond to maxSel selected covariates
      ind <- sapply(maxSel,function(x){which(ResElnet[[i]][[j]]$df==x)[1]})
      lambdavalues <- rep(NaN,length(maxSel))
      for(msel in maxSel){
        fsel <- function(lambda){
          fit<-coef.glmnet(ResElnet[[i]][[j]]$fit,s=lambda,exact=T,x=Xstd[subsamples[,j],],y=Y[subsamples[,j]])
          return(sum(fit!=0)-1-msel) #-1 for intercept
        }
        lmax <- max(1,rev(which(ResElnet[[i]][[j]]$fit$df<msel))[1], na.rm=T) #index of lambda where df is just smaller than msel
        lmax <- ResElnet[[i]][[j]]$fit$lambda[lmax]
        lmin <- which(ResElnet[[i]][[j]]$fit$df>msel)[1] #index of lambda where df is just larger than msel
        lmin <- ResElnet[[i]][[j]]$fit$lambda[lmin]
        if(is.na(lmin)) lmin <- min(ResElnet[[i]][[j]]$fit$lambda)/10^3
        if(is.na(lmax)) lmax <- max(ResElnet[[i]][[j]]$fit$lambda)*10^3
        
        lambdavalues[which(maxSel==msel)] <- uniroot(fsel,interval = c(lmin/2, lmax*2),maxiter = 400,tol=10^(-10))$root
      }
      ResElnet[[i]][[j]]$betaEN<-coef.glmnet(ResElnet[[i]][[j]]$fit,s=lambdavalues,exact=T,x=Xstd[subsamples[,j],],y=Y[subsamples[,j]])
      ResElnet[[i]][[j]]$a0 <- ResElnet[[i]][[j]]$betaEN[1,]
      ResElnet[[i]][[j]]$betaEN <- ResElnet[[i]][[j]]$betaEN[-1,]
      ResElnet[[i]][[j]]$time <-  proc.time()[[3]]-tic
      
      which25EN[[i]][[j]]<-which(ResElnet[[i]][[j]]$betaEN!=0)
    }
    
    write(paste(Sys.time(),"Subsample ",j," done, elastic net",sep=""),file=logname,append=T)
  }
  
  names(which25EN)<-paste("Elnet.alpha.",alp,sep="")
  
  which25 <- c(which25ecpc,which25EN)
  
  save(Res,ResElnet,which25,dfGrps,file = "SubsamplingVerlaat")
}

#Do posterior selection for ecpc on subsamples----
postselection <- "elnet+dense"
maxSel2 <- 50
load(paste(pathResults,"SubsamplingVerlaat",sep=""))

if(0){
  load(paste(pathResults,"which25ecpc",sep=""))
  fname <- paste("which",maxSel2,"ecpcA+",postselection,sep="")
  which50ecpc <- list()
  post2 <- list()
  
  for(i in 1:nSmpl){
    post2[[i]] <- postSelect(X=Xstd[subsamples[,i],],Y=Y[subsamples[,i]],
                           beta=Res[[i]]$beta,intrcpt=Res[[i]]$intercept,
                           postselection=postselection,maxsel=maxSel2, 
                           penalties=Res[[i]]$penalties,model="logistic",
                           tauglobal=Res[[i]]$tauglobal,
                           X2=Xstd[-subsamples[,i],],Y2=Y[-subsamples[,i]])
    
    
    which50ecpc$ecpc[[i]] <- apply(post2[[i]]$betaPost,2,function(x){which(x!=0)})
    
  }
  save(which50ecpc,post2,file=fname) #this comparison
  
  which50EN <-list()
  ResElnet <- list()
  for(i in 1:length(alp)){
    which50EN[[i]]<-list()
    ResElnet[[i]] <- list()
    for(j in 1:nSmpl){
      ResElnet[[i]][[j]]<-list()
      
      tic<-proc.time()[[3]]
      ResElnet[[i]][[j]]$fit<-glmnet(Xstd[subsamples[,j],],Y[subsamples[,j]],alpha=alp[i],family="binomial",
                                     standardize = F,intercept=T) 
      
      #then find lambda values that correspond to maxSel selected covariates
      ind <- sapply(maxSel2,function(x){which(ResElnet[[i]][[j]]$df==x)[1]})
      lambdavalues <- rep(NaN,length(maxSel2))
      for(msel in maxSel2){
        fsel <- function(lambda){
          fit<-coef.glmnet(ResElnet[[i]][[j]]$fit,s=lambda,exact=T,x=Xstd[subsamples[,j],],y=Y[subsamples[,j]])
          return(sum(fit!=0)-1-msel) #-1 for intercept
        }
        lmax <- max(1,rev(which(ResElnet[[i]][[j]]$fit$df<msel))[1], na.rm=T) #index of lambda where df is just smaller than msel
        lmax <- ResElnet[[i]][[j]]$fit$lambda[lmax]
        lmin <- which(ResElnet[[i]][[j]]$fit$df>msel)[1] #index of lambda where df is just larger than msel
        lmin <- ResElnet[[i]][[j]]$fit$lambda[lmin]
        if(is.na(lmin)) lmin <- min(ResElnet[[i]][[j]]$fit$lambda)/10^3
        if(is.na(lmax)) lmax <- max(ResElnet[[i]][[j]]$fit$lambda)*10^3
        
        lambdavalues[which(maxSel2==msel)] <- uniroot(fsel,interval = c(lmin/2, lmax*2),maxiter = 400,tol=10^(-10))$root
      }
      ResElnet[[i]][[j]]$betaEN<-coef.glmnet(ResElnet[[i]][[j]]$fit,s=lambdavalues,exact=T,x=Xstd[subsamples[,j],],y=Y[subsamples[,j]])
      ResElnet[[i]][[j]]$a0 <- ResElnet[[i]][[j]]$betaEN[1,]
      ResElnet[[i]][[j]]$betaEN <- ResElnet[[i]][[j]]$betaEN[-1,]
      ResElnet[[i]][[j]]$time <-  proc.time()[[3]]-tic
      
      which50EN[[i]][[j]]<-which(ResElnet[[i]][[j]]$betaEN!=0)
    }
  }
  
  names(which50EN)<-paste("Elnet.alpha.",alp,sep="")
  
  which50 <- c(which50ecpc,which50EN)
  
  save(which50,post2,ResElnet,file = "SubsamplingVerlaatPosthoc50")
}




#Do CV for ecpc with sds as third grouping ----
fname <- "CVMethylation_ecpc+sds" #file in which the results are stored

if(0){
  grpsno <- c(unlist(sapply(GroupingsAll,function(x){1:length(x)}))) #vector with group numbers in all groupings
  grpngsno <- c(unlist(sapply(1:length(GroupingsAll),function(i){rep(i,length(GroupingsAll[[i]]))}))) #vector with grouping numbers
  dfGrps<-data.frame() #data frame in which group and grouping weights are stored
  dfBeta <- data.frame()
  df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
  Res<-list() #list in which raw output of ecpc is stored (e.g. estimated regression coefficients)
  #for(i in 1:nfolds){
  finalMatrix <- foreach(i=1:nfolds, .combine=rbind,
                         .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                       "CVXR","GRridge","expm","Rsolnp","ecpc")) %dopar% {
       tic<-proc.time()[[3]]
       Res[[i]]<-ecpc(Y[-folds2[[i]]],Xstd[-folds2[[i]],],GroupingsAll,hypershrinkage=hypershrinkage,
                      Y2=Y[folds2[[i]]],X2=Xstd[folds2[[i]],],lambda=lambdas[i],
                      groupings.grouplvl=groupings.grouplvl,postselection="elnet+dense",maxsel=maxSel)
       Res[[i]]$timeGR <- proc.time()[[3]]-tic
       df2<-data.frame("Ypred"=c(Res[[i]]$Ypred,Res[[i]]$Ypredridge,c(Res[[i]]$YpredPost)))
       df2$Method <- rep(c("ecpc","ordinary.ridge",paste("ecpc",maxSel,"vars",sep="")),each=length(folds2[[i]]))
       df2$NumberSelectedVars <- rep(c(p,p,maxSel),each=length(folds2[[i]]))
       df2$Fold <- i
       df2$Sample <- rep(folds2[[i]],2+length(maxSel))
       df2$Time <-  Res[[i]]$timeGR
       df2$Truth <- rep(Y[folds2[[i]]],2+length(maxSel))
       #df<-rbind(df,df2)
       
       df3<-data.frame("Group"=grpsno,
                       "Grouping"=grpngsno,
                       "Group weight"=Res[[i]]$gamma,
                       "Grouping weight"=Res[[i]]$w[grpngsno])
       df3$Tau.ecpc <- Res[[i]]$tauglobal #global tau
       df3$Tau.ridge <- 1/Res[[i]]$lambdaridge #ordinary ridge tau
       df3$Method <- "ecpc"
       df3$Fold <- i
       #dfGrps<-rbind(dfGrps,df3)
       
       whichPost2 <- apply(Res[[i]]$betasPost,2,function(x){which(x!=0)})
       betasPost2 <- unlist(lapply(1:dim(Res[[i]]$betasPost)[2],function(x){Res[[i]]$betasPost[whichPost2[[x]],x]}))
       setting2 <- lapply(1:length(whichPost2),function(x){rep(x,length(whichPost2[[x]]))})
       MaxSel2 <- lapply(1:length(whichPost2),function(x){rep(length(whichPost2[[x]]),length(whichPost2[[x]]))})
       dfBeta2<-data.frame("Betas"=c(betasPost2),
                           "whichPost"=c(unlist(whichPost2)),
                           "Setting"=c(unlist(setting2)),
                           "MaxSel"=c(unlist(MaxSel2)),
                           "Method"=c(rep("ecpc+selection",length(betasPost2))))
       
       dfBeta2$Fold <- i
       #dfBeta<-rbind(dfBeta,dfBeta2)
       
       write(paste(Sys.time(),"fold",i,"of",nfolds,"done"),file=logname,append=T)
       
       list("Res"=Res,"df"=df2,"dfGrps"=df3,"dfBeta"=dfBeta2)
     }
  
  Res <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]][[i]])
  df2 <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]])
  dfGrps2 <- lapply(1:nfolds,function(i) finalMatrix[i,3][[1]])
  dfBeta2 <- lapply(1:nfolds,function(i) finalMatrix[i,4][[1]])
  df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
  dfGrps <- dfGrps2[[1]]; for(i in 2:nfolds) dfGrps <- rbind(dfGrps,dfGrps2[[i]])
  dfBeta <- dfBeta2[[1]]; for(i in 2:nfolds) dfBeta <- rbind(dfBeta,dfBeta2[[i]])
  save(Res,df,dfGrps,dfBeta,file="CVVerlaat2")
  stopCluster(cl)
  
  Summdf <- df %>% group_by(Method) %>% summarise(Brier = mean((Ypred-Truth)^2),
                                                  CVLL = sum(log(Ypred)*Truth+log(1-Ypred)*(1-Truth)),
                                                  NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  df$Method<-as.factor(df$Method)
  dfGrps$Group<-as.factor(dfGrps$Group)
  
  dfROC<-data.frame()
  for(i in levels(df$Method)){
    temp<-data.frame()
    cutoffs<-rev(seq(0,1,by=0.001))
    rocGR <- roc(probs=df$Ypred[df$Method==i],true=df$Truth[df$Method==i],cutoffs=cutoffs)
    temp<-data.frame("FPR"=rocGR[1,],"TPR"=rocGR[2,],"Accuracy"=rocGR[3,])
    temp$Method <- i
    temp$AUC<-c(auc(rocGR))
    temp$NumberSelectedVars<-mean(df$NumberSelectedVars[df$Method==i])
    dfROC<-rbind(dfROC,temp)
  }
  dfAUC <- dfROC %>% group_by(Method) %>% summarise(AUC=mean(AUC),
                                                    NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  
  #plot(dfAUC$NumberSelectedVars,dfAUC$AUC)
  
  if(all(Summdf$NumberSelectedVars==dfAUC$NumberSelectedVars)){ #add AUC to the summary stastistics data frame
    Summdf$AUC <- dfAUC$AUC
  }
  dfGrps$Group<-as.factor(dfGrps$Group)
  save(df,dfGrps,Summdf,dfROC,dfAUC,Res,dfBeta,file=fname) #this comparison 
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
display.brewer.all(3,colorblindFriendly=T)
colpal <- "Dark2"
colsfill <- brewer.pal(3,"Dark2")[1:2]
colsAUC <- brewer.pal(8,"Dark2")

#Load results, given in files:----
load(paste(pathResults,"CVMethylation_ecpcA",sep="")) #df, dfGrps, dfBeta 
dfAUC$Method[!(dfAUC$Method%in%c("ecpc","ordinary.ridge"))] <- "ecpc+selection" #"ecpc+dense+selection"
dfAUC$Sparsity <- "dense"
dfAUC$Sparsity[dfAUC$Method=="ecpc+selection"] <- "sparse"
dfAUC$Codata <- "All"
dfAUC$Codata[dfAUC$Method=="ordinary.ridge"] <- "NA"
dfAUC2<- dfAUC

load(paste(pathResults,"CVMethylation_ecpcB",sep="")) #df, dfGrps, dfBeta
nrSelectedVars <- sapply(Res,function(x){sum(x$beta!=0)})
dfAUC3 <- df[df$Method=="ecpc",] %>% distinct(Fold, .keep_all=T) %>% select(NumberSelectedVars,Method,Fold) %>%
  mutate(AUC=dfAUC$AUC[dfAUC$Method=="ecpc"],Sparsity="group sparse",Codata="All") %>% ungroup()
dfAUC3$NumberSelectedVars <- nrSelectedVars[dfAUC3$Fold]
dfAUC3 <- dfAUC3 %>% select(-Fold)
dfAUC3$Method <- as.factor("ecpc")
dfAUC2 <- rbind(dfAUC2,dfAUC3)

load(paste(pathResults,"CVMethylation_ElNet",sep="")) #dfElNet, dfBetaElNet
dfAUCElNet$Method <- paste("ElNet.","alpha.",dfAUCElNet$Alpha,sep="")
dfAUCElNet$Sparsity <- "sparse"
dfAUCElNet$Codata <- "NA"
dfAUCElNet <- dfAUCElNet %>% select(-c("Alpha","Setting"))
dfAUC2 <- rbind(dfAUC2,dfAUCElNet)

load(paste(pathResults,"CVMethylation_GRridge",sep="")) #dfGRridge, dfGrpsGRridge, dfBetaGRridge
dfAUCGRridge$Method[!(dfAUCGRridge$Method=="GRridge")] <- "GRridge+selection"
dfAUCGRridge$Sparsity <- "dense"
dfAUCGRridge$Sparsity[dfAUCGRridge$Method=="GRridge+selection"] <- "sparse"
dfAUCGRridge$Codata <- "All"
dfAUC2<- rbind(dfAUC2,dfAUCGRridge)

load(paste(pathResults,"CVMethylation_GL",sep="")) #dfGL, dfGrpsGL
dfAUCGL$Sparsity <- "group sparse"
dfAUCGL$Codata <- "All"
dfAUC2 <- rbind(dfAUC2,dfAUCGL)

load(paste(pathResults,"CVMethylation_HierLasso",sep="")) #dfHL, dfGrpsHL
dfAUCHL$Sparsity <- "group sparse"
dfAUCHL$Codata <- "All"
dfAUC2 <- rbind(dfAUC2,dfAUCHL)

load(paste(pathResults,"CVMethylation_RF",sep="")) #dfRF
dfAUCRF$Sparsity <- "dense"
dfAUCRF$Codata <- "NA"
dfAUC2 <- rbind(dfAUC2,dfAUCRF)

dfAUC2$Sparsity <- factor(dfAUC2$Sparsity,levels=c("dense","group sparse","sparse"))
dfAUC2$Codata <- as.factor(dfAUC2$Codata)
dfAUC2$Method <- factor(dfAUC2$Method,levels = c("ecpc","ecpc+selection","GRridge","GRridge+selection",
                                                 "ordinary.ridge","ElNet.alpha.0.3",
                                                 "ElNet.alpha.0.8","GroupLasso",
                                                 "HierLasso","RandomForest"))
lvlsMethod <- c("ecpc","ecpc+selection","GRridge","GRridge+selection",
          "ordinary ridge",bquote("elastic net "*alpha==.3),bquote("elastic net "*alpha==.8),"group lasso",
          "hierarchical lasso","random forest")

#Plot AUC for different methods, dense----
lims <- range(dfAUC2$AUC)
labels<-which(levels(dfAUC2$Method)%in%unique(dfAUC2$Method[dfAUC2$Sparsity=="dense"]))

figname<-paste("FigVerlaatAUCDense.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
# png(filename = figname,
#     width=wdth,
#     height=hght)
#barplot
ggplot(dfAUC2[dfAUC2$Sparsity=="dense",])+
  geom_col(aes(x=Method,y=AUC,fill=Method))+
  scale_fill_manual(values = colsAUC[1:4],labels=lvlsMethod[labels])+
  facet_grid(.~Sparsity,scales="free_x")+
  coord_cartesian(ylim=lims)+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
#lineplot
# ggplot(dfAUC2[dfAUC2$Sparsity=="dense",])+
#   geom_hline(aes(yintercept=AUC,col=Method),size=ls,alpha=0.6)+
#   geom_point(aes(x=NumberSelectedVars,y=AUC,col=Method),size=ps,stroke=strk)+
#   scale_x_continuous(breaks=p)+
#   scale_color_manual(values = colsAUC[1:4])+
#   facet_grid(.~Sparsity,scales="free_x")+
#   labs(x="# parameters")+
#   ylim(lims)+
#   #xlim(c(0,100))+
#   theme_bw()+
#   theme(axis.text.x=element_text(size=ts),
#         axis.text.y=element_text(size=ts),
#         axis.title.x=element_text(size=ts+2),
#         axis.title.y=element_text(size=ts+2),
#         legend.text=element_text(size=ts),
#         legend.title=element_text(size=ts+2),
#         strip.text=element_text(size=ts))#,
#legend.position="bottom")+
#guides(col=guide_legend(nrow=3,byrow=F))
dev.off()

#Plot AUC for different methods, group sparse----
lims <- range(dfAUC2$AUC)
labels<-which(levels(dfAUC2$Method)%in%unique(dfAUC2$Method[dfAUC2$Sparsity=="group sparse"]))

figname<-paste("FigVerlaatAUCGroupSparse.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
# png(filename = figname,
#     width=wdth,
#     height=hght)
ggplot(dfAUC2[dfAUC2$Sparsity=="group sparse",])+
  aes(y=NumberSelectedVars,x=AUC,col=Method)+
  geom_boxplot(position="identity",width=0.01,size=ls/2,alpha=0.2)+
  geom_point(alpha=1,size=ps)+
  coord_flip()+
  scale_color_manual(values = colsAUC[c(1,7,8)],labels=lvlsMethod[labels])+
  facet_grid(.~Sparsity,scales="free_x")+
  labs(y="# selected covariates")+
  #xlim(c(0,100))+
  xlim(lims)+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
#legend.position="bottom")+
#guides(shape=guide_legend(nrow=2,byrow=F),col=guide_legend(nrow=3,byrow=F))
dev.off()

#Plot AUC for different methods, covariate sparse----
lims <- range(dfAUC2$AUC)
labels<-which(levels(dfAUC2$Method)%in%unique(dfAUC2$Method[dfAUC2$Sparsity=="sparse"]))

figname<-paste("FigVerlaatAUCCovSparse.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
# png(filename = figname,
#     width=wdth,
#     height=hght)
ggplot(dfAUC2[dfAUC2$Sparsity=="sparse",])+
  aes(x=NumberSelectedVars,y=AUC,col=Method)+
  geom_line(data=dfAUC2[dfAUC2$Sparsity=="sparse",],
            aes(col=Method),size=ls,alpha=0.8)+ #ecpc+selection, elastic net and GRridge+selection
  geom_point(size=ps/2,stroke=strk/2)+
  facet_grid(.~Sparsity,scales="free_x")+
  #geom_hline(data=dfAUCDense,aes(yintercept=AUC,col=Method,linetype=Sparsity),size=ls)+ #ecpc, GRridge, HL
  #scale_linetype_manual(name="Sparsity",values=c(2,1))+
  #scale_color_discrete(name="Shrinkage method")+
  scale_color_manual(values = colsAUC[c(1,2,5,6)],labels=lvlsMethod[labels])+
  labs(x="# selected covariates")+
  xlim(c(0,100))+
  #scale_x_continuous(trans='log10')+
  ylim(lims)+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
#legend.position="bottom")+
#guides(col=guide_legend(nrow=3,byrow=F))
dev.off()

#Plot AUC for ecpc with different posterior selection methods----
load(paste(pathResults,"CVMethylation_ecpcA",sep="")) #df, dfGrps, dfBeta 
postselection <- c("BRmarginal,dense","BRmarginal,sparse","elnet,dense","elnet,sparse","DSS")
dfAUCPost2 <- data.frame()
for(sel in postselection){
  fnamePost <- paste("CVMethylation_ecpcA+",sel,sep="")
  load(paste(pathResults,fnamePost,sep=""))
  dfAUCPost$Method2 <- factor(paste("ecpc+",sel,sep=""))
  dfAUCPost2 <- rbind(dfAUCPost2,dfAUCPost)
}
dfAUCPost<-dfAUCPost2

#rename factor levels method
lvlsPost <- levels(dfAUCPost$Method2)
dfAUCPost$Method2 <- factor(dfAUCPost$Method2,levels=levels(dfAUCPost$Method2),labels=lvlsPost)
labels<-which(levels(dfAUC2$Method)%in%unique(dfAUC2$Method[dfAUC2$Sparsity=="dense"]))

figname<-paste("FigVerlaatAUCvsNumberParametersPost.pdf",sep="")
pdf(width = wdthpdf*1.25, height = hghtpdf,
    file = figname)
# png(file=figname,
#     width=wdth*1.25,height=hght)
ggplot(dfAUC2)+
  aes(x=NumberSelectedVars,y=AUC,col=Method)+
  #geom_line(aes(col=Method,linetype=Sparsity),size=ls)+ #ecpc+selection, elastic net and GRridge+selection
  geom_hline(data=dfAUC2[dfAUC2$Sparsity=="dense",],
             aes(yintercept=AUC,linetype=Method),col="black",size=ls/1.5,alpha=0.5)+ #ecpc, GRridge, HL
  geom_line(data=dfAUCPost,aes(x=NumberSelectedVars,y=AUC,col=Method2),size=ls,alpha=0.2)+ #ecpc+other post hoc
  geom_point(data=dfAUCPost,aes(x=NumberSelectedVars,y=AUC,col=Method2),size=sz,stroke=strk)+ 
  #scale_linetype_manual(name="Sparsity",values=c(2,1))+
  #scale_color_discrete(name="Post-hoc selection method")+
  scale_color_brewer(name="Posterior selection method",palette=palette)+
  scale_linetype_discrete(name="Dense method",labels=lvlsMethod[labels])+
  labs(x="# selected covariates")+
  #xlim(c(0,100))+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))
dev.off()

#Plot visualisation of a grouping in a graph----
labelsz <- 5.4
ls<-1
load(paste(pathResults,"CVMethylation_ecpcA",sep="")) #df, dfGrps, dfBeta, Res
GroupingsAll <- list(CpGislands=GroupingCpG,
                     Pvalues = GroupingPvals)
groupings.grouplvl <- list(NULL,HierarchyPvals)
grpngsno <- c(unlist(sapply(1:length(GroupingsAll),function(i){rep(i,length(GroupingsAll[[i]]))})))
  
i <- 2  #grouping number
j <- 1 #fold
PlotGraphGrouping <- visualiseGrouping(GroupingsAll[[i]],Res[[j]]$gamma[grpngsno==i],grouping.grouplvl = groupings.grouplvl[[i]],
                                       nodeSize=labelsz,ls=ls)

figname<-paste("FigVerlaatGrouping",i,"withweights.pdf",sep="")
pdf(width = wdthpdf/1.4, height = hghtpdf/1.4,
    file = figname)
# png(file=figname,
#     width=wdth/1.4,height=hght/1.4)
PlotGraphGrouping+
  theme(legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts),
        legend.position = "none")
dev.off()

PlotGraphGrouping2 <- visualiseGrouping(GroupingsAll[[i]],grouping.grouplvl = groupings.grouplvl[[i]],nodeSize=labelsz,ls=ls)
figname<-paste("FigVerlaatGrouping",i,".pdf",sep="")
pdf(width = wdthpdf/1.4, height = hghtpdf/1.4/2,
    file = figname)
# png(file=figname,
#     width=wdth/1.4,height=hght/1.4)
PlotGraphGrouping2+
  theme(legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))
dev.off()

#Plot cross-validated group weight estimates----
library(scales)
load(paste(pathResults,"CVMethylation_ecpcA",sep="")) #df, dfGrps, dfBeta, Res
GroupingsAll <- list(CpGislands=GroupingCpG,
                     Pvalues = GroupingPvals)
groupings.grouplvl <- list(NULL,HierarchyPvals)
values <- list(NULL,pvalFarkas)

names(GroupingsAll) <- c("CpG-islands","p-values")

i <- 2 #grouping number
widthBoxplot <- 0.3; ps<-4
if(i==2) widthBoxplot <- 0.5; ps<-3
dfGrps2 <- dfGrps[dfGrps$Method=="ecpc"&dfGrps$Grouping==i,] %>% select(c(Group.weight,Group,Fold))
PlotGroupweights <- visualiseGroupweights(dfGrps=dfGrps2,Grouping=GroupingsAll[[i]],grouping.grouplvl = groupings.grouplvl[[i]],values=values[[i]],
                                          boxplot=F,ps=ps,ls=ls,widthBoxplot = widthBoxplot)

nCol <- length(levels(PlotGroupweights$data[[PlotGroupweights$labels$colour]])) #number of groups/colors
colsGroupings <- brewer.pal(max(3,max(dfGrps$Grouping)),"Dark2")[1:max(dfGrps$Grouping)]
limsCol <- seq_gradient_pal(low="white",high=colsGroupings[i])((0:2)/2)[2]
limsCol[2]<-seq_gradient_pal(low=colsGroupings[i],high="black")((0:2)/2)[2]
#limsCol<-seq_gradient_pal(low=colsGroupings[i],high="black")((0:2)/3)[1:2]

lims <- range(c(sapply(Res,function(x) range(x$gamma))))
lims<-c(0,51)
figname<-paste("FigVerlaatGroupweightsGrouping",i,".pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
# png(file=figname,
#     width=wdth,height=hght)
plt<-PlotGroupweights+
  theme_bw()+
  scale_color_manual(values=seq_gradient_pal(low=limsCol[1],high="black")((0:(nCol-1))/(nCol)),name="Group",guide=F)+
  labs(title=paste("Grouping",names(GroupingsAll)[i]))+
  ylim(ylim=lims)+
  theme(title=element_text(size=ts+2),
        axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
#legend.position="bottom")+
#guides(col=guide_legend(nrow=3,byrow=F))
if(i==2){
  plt <- plt+scale_color_brewer(name="Group",palette="Dark2")+
    scale_x_continuous(trans='log10')+
    labs(x="p-values")
}
plt
dev.off()

#Plot cross-validated groupings weight estimates----
load(paste(pathResults,"CVMethylation_ecpcA",sep="")) #df, dfGrps, dfBeta, Res
GroupingsAll <- list(CpGislands=GroupingCpG,
                     Pvalues = GroupingPvals)
dfGrps3 <- dfGrps[dfGrps$Group==1,]
PlotGroupingWeights <- visualiseGroupingweights(dfGrps=dfGrps3,GroupingNames=names(GroupingsAll),
                                                boxplot=F,ps=5,width=0.25)

figname<-paste("FigVerlaatGroupingweights.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
# png(file=figname,
#     width=wdth,height=hght)
PlotGroupingWeights+
  scale_color_brewer(name="Grouping",palette="Dark2")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))
dev.off()



#Plot histogram of beta estimates ridge and ecpc----
load(paste(pathResults,"CVMethylation_ecpcA",sep="")) #df, dfGrps, dfBeta 

i<-1
dfPlot <- data.frame(Beta=c(Res[[i]]$beta,Res[[i]]$betaridge),
                     Method=as.factor(rep(c("ecpc","ordinary ridge"),each=p)))

q09 <- quantile(abs(dfPlot$Beta),0.9)
ylim <- c(0,12)
figname<-paste("FigVerlaatDensityBetas.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
# png(file=figname,
#     width=wdth,height=hght)
ggplot(dfPlot)+aes(x=abs(Beta))+
  geom_histogram(aes(fill=Method,col=Method,y=..density..),alpha=0.2,position = "identity")+
  #geom_density(aes(col=Method),size=ls)+
  geom_vline(xintercept=q09,linetype=2)+
  xlab(label=parse(text="abs( beta)"))+
  scale_color_brewer(palette=palette)+
  scale_fill_brewer(palette=palette)+
  #xlim(c(q09,NA))+
  #coord_cartesian(xlim=c(q09,NA),ylim=ylim)+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))
dev.off()
#Plot AUC in subsamples----
load("subsamplesVerlaat")
nSmpl<-50

#save(dfGrps,Res,ResElnet,which25,file="SubsamplingVerlaat")
#save(post2,ResElnet,which50,file="SubsamplingVerlaatPosthoc50")
# for(i in 1:nSmpl){
#   post2[[i]]$betaPost<-post2[[i]]$betasEN
#   post2[[i]]$interceptPost <- post2[[i]]$interceptEN
#   ResElnet[[1]][[i]]$betaEN <- ResElnet[[1]][[i]]$betasEN
#   ResElnet[[2]][[i]]$betaEN <- ResElnet[[2]][[i]]$betasEN
# }

if(1){
  figname<-paste("FigVerlaatAUCSubsamples25.pdf",sep="")
  load(paste(pathResults,"SubsamplingVerlaat",sep="")) #25 selected
  brierScore.ecpc <- sapply(1:nSmpl,function(i){
    Ytrue <- Y[-subsamples[,i]] #left out samples
    #X2 <- Xstd[subsamples[,i],] #left out data
    lps <- Xstd[-subsamples[,i],]%*%Res[[i]]$betaPost + Res[[i]]$interceptPost[1]
    Ypred <- 1/(1+exp(-lps))
    return(mean((Ytrue-Ypred)^2))
  })
  auc.ecpc <- sapply(1:nSmpl,function(i){
    Ytrue <- Y[-subsamples[,i]] #left out samples
    #X2 <- Xstd[subsamples[,i],] #left out data
    lps <- Xstd[-subsamples[,i],]%*%Res[[i]]$betaPost + Res[[i]]$interceptPost[1]
    Ypred <- 1/(1+exp(-lps))
    
    cutoffs<-rev(seq(0,1,by=0.001))
    ROC <- roc(probs=Ypred,true=Ytrue,cutoffs=cutoffs)
    return(auc(ROC))
  })
}else{
  figname<-paste("FigVerlaatAUCSubsamples50.pdf",sep="")
  load(paste(pathResults,"SubsamplingVerlaatPosthoc50",sep="")) #50 selected
  brierScore.ecpc <- sapply(1:nSmpl,function(i){
    Ytrue <- Y[-subsamples[,i]] #left out samples
    #X2 <- Xstd[subsamples[,i],] #left out data
    lps <- Xstd[-subsamples[,i],]%*%post2[[i]]$betaPost + post2[[i]]$a0[1]
    Ypred <- 1/(1+exp(-lps))
    return(mean((Ytrue-Ypred)^2))
  })
  auc.ecpc <- sapply(1:nSmpl,function(i){
    Ytrue <- Y[-subsamples[,i]] #left out samples
    #X2 <- Xstd[subsamples[,i],] #left out data
    lps <- Xstd[-subsamples[,i],]%*%post2[[i]]$betaPost + post2[[i]]$a0[1]
    Ypred <- 1/(1+exp(-lps))
    
    cutoffs<-rev(seq(0,1,by=0.001))
    ROC <- roc(probs=Ypred,true=Ytrue,cutoffs=cutoffs)
    return(auc(ROC))
  })
  which25<-which50
}

brierScore.elnet <- lapply(ResElnet,function(ResPost){
  sapply(1:nSmpl,function(i){
    Ytrue <- Y[-subsamples[,i]] #left out samples
    #X2 <- Xstd[-subsamples[,i],] #left out data
    lps <- Xstd[-subsamples[,i],]%*%ResPost[[i]]$betaEN + ResPost[[i]]$a0
    Ypred <- 1/(1+exp(-lps))
    return(mean((Ytrue-Ypred)^2))
  })
})
auc.elnet <- lapply(ResElnet,function(ResPost){
  sapply(1:nSmpl,function(i){
    Ytrue <- Y[-subsamples[,i]] #left out samples
    #X2 <- Xstd[-subsamples[,i],] #left out data
    lps <- Xstd[-subsamples[,i],]%*%ResPost[[i]]$betaEN + ResPost[[i]]$a0
    Ypred <- 1/(1+exp(-lps))
    
    cutoffs<-rev(seq(0,1,by=0.001))
    ROC <- roc(probs=Ypred,true=Ytrue,cutoffs=cutoffs)
    return(auc(ROC))
  })
})

df <- data.frame("Brierscore"=c(brierScore.ecpc,unlist(brierScore.elnet)),
                 "AUC"=c(auc.ecpc,unlist(auc.elnet)),
                 "Method"=factor(rep(names(which25),each=nSmpl)))
levels(df$Method)
lvlsMethods<- c("ecpc",bquote("elastic net "*alpha==.3),bquote("elastic net "*alpha==.8))

pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
# png(file=figname,
#     width=wdth,height=hght)
p3<-ggplot(df)+
  geom_boxplot(aes(y=AUC,x=Method,fill=Method),alpha=1)+
  scale_fill_brewer(palette=palette,labels=lvlsMethods)+
  #labs(y="Brier score")+
  ylim(c(0.33,0.9))+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))
p3
dev.off()
#Plot overlap in subsamples----
load(paste(pathResults,"SubsamplingVerlaat",sep=""))
sapply(which25[[1]],length)
# dfGrps2 <- dfGrps[[1]]
# for(i in 1:length(dfGrps)) dfGrps2 <- rbind(dfGrps2,dfGrps[[i]])
# dfGrps<-dfGrps2

load(paste(pathResults,"SubsamplingVerlaatPosthoc50",sep="")) #50 selected
which25<-which50

nSmpl<-50
overlap<-lapply(which25,function(whichselected){
  unlist(sapply(1:(nSmpl-1),function(x){
    sapply((x+1):nSmpl,function(y){
      overlap <- sum(whichselected[[x]]%in%whichselected[[y]])
      return(overlap)
    })
  }))
})
comp1 <- unlist(sapply(1:(nSmpl-1),function(x){rep(x,length((x+1):nSmpl))})) #number of bootstrap sample
comp2 <- unlist(sapply(1:(nSmpl-1),function(x){(x+1):nSmpl}))

df <- data.frame("Overlap"=unlist(overlap),"Method"=rep(names(overlap),each=length(overlap[[1]])))
levels(df$Method)
lvlsMethods<- c("ecpc",bquote("elastic net "*alpha==.3),bquote("elastic net "*alpha==.8))
df %>% group_by(Method) %>% summarise("medianOverlap"=median(Overlap))

figname<-paste("FigVerlaatOverlap50.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
# png(file=figname,
#     width=wdth,height=hght)
p3<-ggplot(df)+
  geom_histogram(aes(x=Overlap,fill=Method),alpha=1,position="dodge",center=0,binwidth=1)+
  scale_fill_brewer(palette = palette,labels=lvlsMethods)+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))
p3
dev.off()

#Plot group weights estimates in subsamples----
load(paste(pathResults,"which25ecpc",sep="")) #dfGrps
GroupingsAll <- list(CpGislands=GroupingCpG,
                     Pvalues = GroupingPvals)
groupings.grouplvl <- list(NULL,HierarchyPvals)
values <- list(NULL,pvalFarkas)

i <- 1  #grouping number
dfGrps2 <- dfGrps[dfGrps$Method=="ecpc"&dfGrps$Grouping==i,] %>% select(c(Group.weight,Group,Fold))
PlotGroupweights <- visualiseGroupweights(dfGrps=dfGrps2,Grouping=GroupingsAll[[i]],grouping.grouplvl = groupings.grouplvl[[i]],values=values[[i]])

figname<-paste("FigVerlaatBTGroupweightsGrouping",i,".png",sep="")
png(file=figname,
    width=wdth,height=hght)
PlotGroupweights+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))
dev.off()

#Plot groupings weight estimates in subsamples----
load(paste(pathResults,"which25ecpc",sep="")) #dfGrps
GroupingsAll <- list(CpGislands=GroupingCpG,
                     Pvalues = GroupingPvals)
dfGrps3 <- dfGrps[dfGrps$Group==1,]
PlotGroupingWeights <- visualiseGroupingweights(dfGrps=dfGrps3,GroupingNames=names(GroupingsAll))

figname<-paste("FigVerlaatSubGroupingweights.png",sep="")
png(file=figname,
    width=wdth,height=hght)
PlotGroupingWeights+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))
dev.off()


#Compare AUC ecpc with or without sds----
load(paste(pathResults,"CVMethylation_ecpcA",sep="")) 
dfAUC2<-dfAUC
head(dfAUC2) #ecpc dense: AUC=0.674
dfAUC2[which.max(dfAUC2$AUC),] #max at 4 selected vars: AUC=0.732

load(paste(pathResults,"CVMethylation_ecpc+sds",sep="")) #ecpc with sds as a third grouping
dfAUC2<-dfAUC
head(dfAUC2) #ecpc dense: AUC=0.663
dfAUC2[which.max(dfAUC2$AUC),] #max at 4 selected vars: AUC=0.737

#Compare AUC results graper and ecpc----
load(paste(pathResults,"CVMethylation_ecpcA",sep="")) 
dfAUC2<-dfAUC
head(dfAUC2) #ecpc dense: AUC=0.674
tail(dfAUC2) #ordinary ridge: AUC=0.671
dfAUC2[which.max(dfAUC2$AUC),] #max at 4 selected vars: AUC=0.732

load(paste(pathResults,"CVMethylation_graper",sep="")) #leaves of FDR2 used, sparse
dfAUC #AUC=0.700

load(paste(pathResults,"CVMethylation_graperA",sep="")) #leaves of FDR2 used, dense
dfAUC #AUC=709

