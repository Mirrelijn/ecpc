###Simulation study
#compare ecpc, ordinary ridge and GRridge
#Co-data: random or informative
#Goal: show that ecpc with hypershrinkage shrinks group weights for random co-data
#      and shrinks little for informative co-data, improving predictions
#NOTE 1: in Rstudio; press ALT+o for overview of code section headings
#NOTE 2: if whole file is run, analysis is skipped,
#      but results are loaded and plots are saved to the working directory

#Load libraries----
setwd("...") #directory in which results and figures are saved
library(ecpc)
# load libraries required for ecpc
library(MASS)
library(penalized)
library(glmnet)
library(mvtnorm)
library(gglasso)
library(mgcv)
library(CVXR)
library(Rsolnp)
library(invgamma)
library(foreach)
library(doParallel)
library(Matrix)
library(GRridge)

#load libraries needed for storing and plotting results
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

#optional: setup parallel backend to use many processors
# cores=detectCores()
if(0){
  cl <- makeCluster(5) #not to overload your computer
  registerDoParallel(cl)
}

#Set simulation variables ----
nSim<-100 #number of simulated test and training data sets
n<-100 #sample size of each training data set
n2<-100 #sample size of each test data set
p<-300 #number of covariates
#simulate betas i.i.d. from beta_k~N(mean=0,sd=sqrt(0.1)), i.e. use one group
muGrp1<-c(0) #group mean
varGrp1<-c(0.1) #group variance
indT1<-rep(1,p) #vector with group numbers; all p beta in group 1

#Generate/load data ----
set.seed(1027439)
Dat1<-list()
for(i in 1:nSim){
  Dat1[[i]]<-simDat(n,p,n2,muGrp1,varGrp1,indT1,sigma=1,model='linear',flag=F) #simulate 100 test and train data sets and betas
}
save(Dat1,muGrp1,varGrp1,indT1,file="SimData1")
#load("SimData1")

#Define range of number of groups used in groupings ----
minGroupSize <- 10 #set minimum group size
maxNoGrps <- floor(p/minGroupSize) #compute corresponding maximum number of groups
temp <- rep(1:9,ceiling(log(maxNoGrps)/log(10))+1)*(10^rep(0:ceiling(log(maxNoGrps)/log(10)),each=9))
G <- (1:maxNoGrps)[temp[temp<=maxNoGrps]] #take number of groups increasing in logarithmic scale (as we want large enough differences in resulting group sizes:)
print("Number of groups tested in analysis:"); G; print("Corresponding average group size:"); p/G

#Make random grouping ----
set.seed(74743)
rankRandom <- sample(1:p,p,replace = F) #ranking used for random groups
ind <- list() #first make list with for each number of groups, a px1 vector with group number
for(g in 1:length(G)){
  size <- floor(p/G[g]) #group size
  rest <- p- size*G[g] #first rest number of groups obtain 1 covariate more than size
  if(rest>0){
    ind[[g]]<-c(rep(1:rest,each=(size+1)),rep((rest+1):G[g],each=size))
  }else{
    ind[[g]]<-rep(1:G[g],each=size)
  }
}

#for each number of groups G[g], create list element with list of G[g] random groups
#e.g. indRandom[[3]][[1]] is the vector with the covariate indices belonging to group 1 in the grouping with G[3] random groups
indRandom <- lapply(ind,function(grps){lapply(1:max(grps),function(x){rankRandom[which(grps==x)]})})

#Keep track of progress in log file----
logname <- "logSimStudy.txt" #
write(paste(Sys.time(),"Log simulation study"),file=logname,append=T)


#Fit Ecpc & ordinary ridge: train and test models on each of the nSim simulated data sets ====
hypershrinkageAll <- c("none","ridge")
hypershrinkage <- hypershrinkageAll[2]
dfG<-data.frame() #create empty data frame to store results for various number of covariate groups

if(0){
  #i<-5
  #for(i in 1:nSim){
  finalMatrix <- foreach(i=1:nSim, .combine=rbind,
                         .packages = c("glmnet","penalized","mvtnorm","gglasso","mgcv",
                                       "CVXR","GRridge","expm","Rsolnp","ecpc")) %dopar% {
       Y<-Dat1[[i]]$Y
       Y2<-Dat1[[i]]$Y2
       X<-Dat1[[i]]$Xctd
       X2<-Dat1[[i]]$X2ctd
       rankBeta<-order(abs(Dat1[[i]]$beta)) #betas ranked on increasing magnitude
       
       dfG2 <- data.frame()
       resRank<-list()
       resRandom<-list()
       ###fit for all groups, for random and informative co-data ----
       for(g in 1:length(G)){
         indRank <- lapply(1:max(ind[[g]]),function(x){rankBeta[which(ind[[g]]==x)]}) #list with G[g] rank-based groups
         
         #No extra shrinkage, or inverse gamma with constraints
         #rank-based groups
         write(paste(Sys.time(),"Rank-based groups:", hypershrinkage),file=logname,append=T)
         resRank[[g]]<-ecpc(Y,X,groupings=list(indRank),Y2=Y2,X2=X2,sigmasq=1,hypershrinkage=hypershrinkage) #informative grouping
         
         #random groups
         write(paste(Sys.time(),"Random groups:", hypershrinkage),file=logname,append=T)
         resRandom[[g]]<-ecpc(Y,X,groupings=list(indRandom[[g]]),Y2=Y2,X2=X2,sigmasq=1,hypershrinkage=hypershrinkage) #random grouping
         
         temp <- data.frame("MSE"=c(resRank[[g]]$MSEecpc[2],resRandom[[g]]$MSEecpc[2]),
                            "Grouping"=rep(c("Rank-based","Random"),each=1),
                            "Method"=rep(c(hypershrinkage),2))
         temp$G <- G[g]
         temp$Sim <- i
         dfG2 <- rbind(dfG2,temp)
         
         write(paste(Sys.time(),"Simulated data set",i,"of",nSim,", number of groups",G[g],"done"),file=logname,append=T)
         
       }
       temp <- data.frame("MSE"=c(resRank[[g]]$MSEridge,resRandom[[g]]$MSEridge),
                          "Grouping"=c("Rank-based","Random"),
                          "Method"=rep("ordinary ridge",2))
       temp$G <- G[g]
       temp$Sim <- i
       dfG2 <- rbind(dfG2,temp)
       
       list("dfG2"=dfG2,"resRank"=resRank,"resRandom"=resRandom)
     }
  
  dfG2 <- lapply(1:nSim,function(i) finalMatrix[i,1][[1]])
  resRank <- lapply(1:nSim,function(i) finalMatrix[i,2][[1]])
  resRandom <- lapply(1:nSim,function(i) finalMatrix[i,3][[1]])
  dfG <- dfG2[[1]]; for(i in 2:nSim) dfG <- rbind(dfG,dfG2[[i]])
  fname <- paste("ResSimDat1_ecpc_",hypershrinkage,sep="")
  save(finalMatrix,resRank,resRandom,dfG,file=fname)
  stopCluster(cl)
  
  
  #Fit GRridge: train and test models on each of the nSim simulated data sets ====
  for(i in 1:nSim){
    Y<-Dat1[[i]]$Y
    Y2<-Dat1[[i]]$Y2
    X<-Dat1[[i]]$Xctd
    X2<-Dat1[[i]]$X2ctd
    rankBeta<-order(abs(Dat1[[i]]$beta)) #betas ranked on increasing magnitude
    
    ###fit for all groups, for random and informative co-data ----
    for(g in 1:length(G)){
      indRank <- lapply(1:max(ind[[g]]),function(x){rankBeta[which(ind[[g]]==x)]}) #list with G[g] rank-based groups
      
      #No extra shrinkage, or inverse gamma with constraints
      #rank-based groups
      write(paste(Sys.time(),"Rank-based groups: GRridge"),file=logname,append=T)
      resRankGR<-grridge(t(X),c(Y),partitions=list(indRank),standardizeX=F)
      Ypred<-X2%*%resRankGR$betas+resRankGR$predobj$GroupRegul$a0
      MSErankGR <- mean((Y2-Ypred)^2)
      
      #random groups
      write(paste(Sys.time(),"Random groups: GRridge"),file=logname,append=T)
      resRandomGR<-grridge(t(X),c(Y),partitions=list(indRandom[[g]]))
      Ypred<-X2%*%resRandomGR$betas+resRandomGR$predobj$GroupRegul$a0
      MSErandomGR <- mean((Y2-Ypred)^2)
      
      temp <- data.frame("MSE"=c(MSErankGR,MSErandomGR),
                         "Grouping"=rep(c("Rank-based","Random"),each=1),
                         "Method"=rep(c("GRridge"),2))
      temp$G <- G[g]
      temp$Sim <- i
      dfG <- rbind(dfG,temp)
      
      write(paste(Sys.time(),"GRridge, simulated data set",i,"of",nSim,", number of groups",G[g],"done"),file=logname,append=T)
      
    }
  }
  
  save(dfG,resRankGR,resRandomGR,file="ResSimDatGRridge") #with other
}

#Plots: general parameters----
wdth<-1200
hght<-450
wdthpdf <- wdth/75
hghtpdf <- hght/75
ts <- 16 #basis text size in figures
ls <- 1.2 #basis line size in figures
ps <- 2 #basis point size in figures
sz <- 2 #point size
strk <- 1.5 #stroke size
# palette <- "Paired"
# cols<-brewer.pal(3,palette)
palette <- "Dark2"
colsAUC <- brewer.pal(8,palette)
colsecpc2 <- seq_gradient_pal(low=colsAUC[1],high="black")((0:3)/3)[2]
cols <- c(colsecpc2,colsAUC)

#Load results----
hypershrinkageAll <- c("none","ridge")

dfG2 <-data.frame()
for(i in 1:length(hypershrinkageAll)){
  fname <- paste("ResSimDat1_ecpc_",hypershrinkageAll[i],sep="")
  load(fname)
  dfG2<-rbind(dfG2,dfG)
}
load("ResSimDatGRridge")
dfG<- rbind(dfG2,dfG)
dfG$Method <- factor(dfG$Method,levels = levels(dfG$Method)[c(1,3,4,2)],
                      labels=c("ecpc without hypershrinkage","ecpc with hypershrinkage","GRridge","ordinary ridge"))


summG <- dfG %>% group_by(Method,Grouping,G) %>% summarise("mean"=mean(MSE),
                                                           "median"=median(MSE),
                                                           "q75"=quantile(MSE,0.75),
                                                           "q25"=quantile(MSE,0.25)) %>% ungroup()
summG$Grouping <- factor(summG$Grouping,levels=levels(summG$Grouping),labels= c("Random","Informative"))
dfG$Grouping <- factor(dfG$Grouping,levels=levels(dfG$Grouping),labels= c("Random","Informative"))

#Plot quantiles of MSE over all nSim simulated data sets ----
# fname<-paste("MSEgroups2.png",sep="")
# png(width = wdth, height = hght,
#     filename = fname)
fname<-paste("MSEgroups2.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = fname)
ggplot(dfG)+aes(x=G)+
  #geom_smooth(aes(color=Method),alpha=0.6)+
  geom_rect(data=summG[summG$Method=="ordinary ridge",],aes(ymin=q25,ymax=q75,fill=Method),
            xmin=min(dfG$G),xmax=max(dfG$G),alpha=0.1)+
  geom_ribbon(data=summG[summG$Method!="ordinary ridge",],aes(ymin=q25,ymax=q75,x=G,fill=Method),alpha=0.1)+
  geom_hline(data=summG[summG$Method=="ordinary ridge",],aes(yintercept=q25,col=Method,linetype=Method),alpha=0.2,size=ls)+ #q25
  geom_hline(data=summG[summG$Method=="ordinary ridge",],aes(yintercept=q75,col=Method,linetype=Method),alpha=0.2,size=ls)+ #q75
  geom_line(data=summG[summG$Method!="ordinary ridge",],aes(x=G,y=q25,color=Method,linetype=Method),alpha=0.2,size=ls)+ #q25
  geom_line(data=summG[summG$Method!="ordinary ridge",],aes(x=G,y=q75,color=Method,linetype=Method),alpha=0.2,size=ls)+ #q75
  geom_hline(data=summG[summG$Method=="ordinary ridge",],aes(yintercept=mean,col=Method,linetype=Method),size=ls*1.2)+ #median
  geom_line(data=summG[summG$Method!="ordinary ridge",],aes(x=G,y=mean,color=Method,linetype=Method),size=ls*1.2)+ #median
  #scale_color_manual(values=c(cols,"black"))+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  scale_linetype_manual(values=c(1,6,2,3))+
  labs(y="MSE")+
  facet_grid(.~Grouping)+
  #facet_grid(Grouping~Method)+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))
dev.off()

#Plot MSE for a couple of simulated data sets --------
i<-c(1:5)
fname<-paste("MSEgroupsNoConstraintsEach.pdf",sep="")
# png(width = wdth*1.2, height = hght,
#     filename = fname)
pdf(width = wdthpdf*1.2, height = hghtpdf,
    file = fname)
ggplot(dfG[dfG$Sim%in%i&dfG$Method!="ordinary ridge",])+aes(x=G,y=MSE,reorder(Sim,G))+
  #geom_smooth(aes(color=Method),alpha=0.6)+
  geom_path(data=dfG[dfG$Sim%in%i&dfG$Method!="ordinary ridge",],aes(group=Sim,col=Method),size=sz,alpha=0.6)+ #mean
  # geom_hline(data=summG[summG$Method=="ordinary ridge",],aes(yintercept=mean),size=sz)+ #mean
  # geom_hline(data=summG[summG$Method=="ordinary ridge",],aes(yintercept=q25),linetype=2,size=sz)+ #q25
  # geom_hline(data=summG[summG$Method=="ordinary ridge",],aes(yintercept=q75),linetype=2,size=sz)+ #q75
  facet_grid(Grouping~Method)+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))
dev.off()
