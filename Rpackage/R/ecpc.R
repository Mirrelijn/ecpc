###ecpc method to learn from co-data:
#Fit a generalised linear (linear, logistic) or Cox model, penalised with adaptive multi-group penalties.
#The method combines empirical Bayes estimation for the group hyperparameters with an extra level of shrinkage
#to be able to handle various co-data, including overlapping groups, hierarchical groups and continuous co-data.

ecpc <- function(Y,X,groupings,groupings.grouplvl=NULL,
                 hypershrinkage,unpen=NULL,intrcpt=T,model,postselection="elnet+dense",maxsel=10,
                 lambda="CV",fold=10,sigmasq=NaN,w=NaN,
                 nsplits=100,weights=T,profplotRSS=F,
                 Y2=NaN,X2=NaN,compare=T,
                 mu=F,normalise=F,silent=F
                 #nIt=1,betaold=NaN
                 ){
  #-1. Description input --------------------------------------------------------------------------
  #
  #Data and co-data:
  # Y: nx1 vector with response data
  # X: nxp matrix with observed data
  # groupings: list of m elements, each element one co-data grouping 
  #            with each grouping a list of groups containing the indices of covariates in that group 
  # groupings.grouplvl: (optional) hierarchical groups define a grouping on group level.
  #            list of m elements (corresponding to groupings), with NULL if there is no structure on group level, or 
  #            with a list of groups containing the indices of groups of covariates in that group
  #
  #Model:
  # hypershrinkage: vector of m strings indicating shrinkage type used in level of extra shrinkage for each of the m groupings, 
  #                either in the the form "type" 
  #                or "type1,type2", in which type1 is used to select groups, and type 2 to estimate selected group parameters
  # unpen: vector with indices of covariates that should not be penalised (TD: adapt mlestlin)
  # intrcpt: T/F to use/not use intercept (always unpenalised)
  # model: type of response Y (linear, logistic or Cox)
  # postselection: T/F if parsimonious model is/is not needed, or string with type of posterior selection method;
  #               "elnet", or "DSS" for corresponding post-selection method used
  # maxsel: vector with maximum number of penalised covariates to be selected (additional to all unpenalised covariates)
  #
  #Global hyperparameter estimation:
  #NOTE: tau^2_{global}=sigma^2/lambda (linear) or 1/lambda (logistic,Cox)
  # lambda: -numeric value for lambda to be used to compute initial beta on which EB estimates are based, or
  #         -"ML" or "CV" if (approximate) ML criterium or cross-validation needs to be used to estimate overall lambda
  # fold: number of folds used in inner CV if tau^_{global}^2 is estimated with CV
  # sigmasq: (linear model) given variance of Y~N(X*beta,sd=sqrt(sigmasq)), else estimated
  # w: (optional) mx1 vector to fix grouping weights at given value.
  #
  #Local hyperparameter estimation:
  # nsplits: number of splits used in RSS criterion in the extra level of shrinkage
  # weights: T/F to use weights in ridge hypershrinkage to correct for group size
  #
  #Optional settings:
  # Y2,X2: (optional) independent data set for performance check 
  # compare: -T if to grridge to be compared with glmnet, with same lambda. 
  #          -if "CV" or "ML", (possibly) different lambda used in glmnet (lambda "ML", "CV" specifies which method is used for grridge lambda)
  #          -for logistic/cox: if "MoM", CV approximation as initial value + MoM for moment iteration on whole group
  # silent: set to T to suppress output messages
  #
  #Experimental settings:
  # mu: T/F to include/exclude group prior means (default F)
  # normalise: T if group variances should be normalised to sum to 1 (default F)
  # nIt: number of times ecpc is iterated (default 1)
  # betaold: if beta known for similar, previous study, can be used as weights for group mean beta_old*mu (default unknown)

  nIt=1;
  betaold=NaN
  #-2. Set-up variables ---------------------------------------------------------------------------
  n <- dim(X)[1] #number of samples
  p <- dim(X)[2] #number of covariates 
  if(!missing(X2)) n2<-dim(X2)[1] #number of samples in independent data set x2 if given
  
  if(missing(model)){
    if(all(is.element(Y,c(0,1))) || is.factor(Y)){
      model <- "logistic" 
    } else if(all(is.numeric(Y)) & !(is.matrix(Y) && dim(Y)[2]==2)){
      model <- "linear"
    }else{
      model <- "cox"
    }
  }
  levelsY<-NaN
  if(is.nan(lambda)) lambda <- "CV"
  if(model=="logistic"){
    levelsY<-cbind(c(0,1),c(0,1))
    if(lambda=="ML"){
      lambda <- "CV"
      if(!silent) print("For logistic model, use CV for overall tau.")
    }
    if(!all(is.element(Y,c(0,1)))){
    oldLevelsY<-levels(Y)
    levels(Y)<-c("0","1")
    Y<-as.numeric(Y)-1
    levelsY<-cbind(oldLevelsY,c(0,1))
    colnames(levelsY)<-c("Old level names","New level names")
    if(!missing(Y2)){
      levels(Y2)<-c("0","1")
      Y2<-as.numeric(Y2)-1
    }
    if(!silent) print("Y is put in 0/1 format, see levelsY in output for new names")
    }
  }
  if(model=='cox'){
    intrcpt <- F
    if(is.nan(tausq)){
      if(is.nan(lambda)) lambda <- "CV"
      if(lambda=="ML"){
        if(!silent) print("For cox model, no ML approximation for overall tau available. Use CV instead.")
        lambda <- "CV"
      }
    }
  }
  switch(model,
         'linear'={
           fml <- 'gaussian'
           sd_y <- sqrt(var(Y)*(n-1)/n)[1]
         },
         'logistic'={
           fml <- 'binomial'
           sd_y <- 1 #do not standardise y in logistic setting
         },
         'cox'={
           fml <- 'cox'
           sd_y <- 1 #do not standardise y in cox regression setting
         }
  )
  mutrgt<-0
  if(mu==F){mu <- 0}else mu <- NaN
  tautrgt<-NaN
  tausq<-NaN
  hyperlambdas<-c(NaN,NaN)
  
  # check whether unpenalised covariates are not in partition
  # and set penalty.factor of unpenalised covariates to 0 for glmnet
  penfctr <- rep(1,p) #factor=1 for penalised covariates
  if(length(unpen)>0){
    penfctr[unpen] <- 0 #factor=0 for unpenalised covariates
    if(any(unlist(groupings)%in%unpen)){
      warning("Unpenalised covariates removed from grouping")
      for(i in 1:length(groupings)){
        for(j in 1:length(groupings[[i]])){
          if(all(groupings[[i]][[j]]%in%unpen)){
            groupings[[i]][[j]] <- NULL #remove whole group if all covariates unpenalised
          }else{
            groupings[[i]][[j]] <- groupings[[i]][[j]][!(groupings[[i]][[j]]%in%unpen)]
          }
        }
      }
    }
  }
  
  #-2.1 Variables describing groups and partition(s) =================================================
  G <- sapply(groupings,length) #1xm vector with G_i, number of groups in partition i
  m <- length(G) #number of partitions
  if(missing(hypershrinkage)){
    hypershrinkage<-rep("ridge",m)
  }
  if(any(grepl("hierLasso",hypershrinkage))){
    if(length(groupings.grouplvl)==0){
      stop("Grouping on group level for hierarchical groups is missing")
    }
    if(!is.list(groupings.grouplvl) | length(groupings.grouplvl)!=m){
      stop("Groupings on group level should be a nested list")
    }
  }
  indGrpsGlobal <- list(1:G[1]) #global group index in case we have multiple partitions
  if(m>1){
    for(i in 2:m){
      indGrpsGlobal[[i]] <- (sum(G[1:(i-1)])+1):sum(G[1:i])
    }
  }
  Kg <- lapply(groupings,function(x)(sapply(x,length))) #m-list with G_i vector of group sizes in partition i
  #ind1<-ind
  
  #ind <- (matrix(1,G,1)%*%ind)==(1:G)#sparse matrix with ij element T if jth element in group i, otherwise F
  i<-unlist(sapply(1:sum(G),function(x){rep(x,unlist(Kg)[x])}))
  j<-unlist(unlist(groupings))
  ind <- sparseMatrix(i,j,x=1) #sparse matrix with ij element 1 if jth element in group i (global index), otherwise 0
  
  Ik <- lapply(1:m,function(i){
    x<-rep(0,sum(G))
    x[(sum(G[1:i-1])+1):sum(G[1:i])]<-1
    as.vector(x%*%ind)}) #list for each partition with px1 vector with number of groups beta_k is in
  #sparse matrix with ij element 1/Ij if beta_j in group i
  
  #make co-data matrix Z (Zt transpose of Z as in paper, with co-data matrices stacked for multiple groupings)
  Zt<-ind; 
  if(G[1]>1){
    Zt[1:G[1],]<-t(t(ind[1:G[1],])/apply(ind[1:G[1],],2,sum))
  }
  if(m>1){
    for(i in 2:m){
      if(G[i]>1){
        Zt[indGrpsGlobal[[i]],]<-t(t(ind[indGrpsGlobal[[i]],])/
                                                     apply(ind[indGrpsGlobal[[i]],],2,sum))
      }
    }
  }
  PenGrps <- as.matrix(Zt%*%t(Zt)) #penalty matrix groups

  #-2.2 Weight variables for extra shrinkage on group parameters =====================================
  # Compute weights and corresponding weight matrix
  #Note: for logistic regression, another, different weight matrix W is defined below
  if(weights){
    weights <- unlist(Kg)
  }else{
    weights <- rep(1,sum(G))
  }
  if(length(weights)==1){Wminhalf<-1/sqrt(weights)}
  else{Wminhalf <- diag(1/sqrt(weights))} #W^{-0.5} (element-wise), with W diagonal matrix with weights for prior parameters

  #-3. The ecpc method possibly using extra shrinkage ---------------------------------------------
  #-3.1 Set-up variables =======================================================================================
  #-3.1.1 Variables adapted for intercept ####################################################################
  # copy X: add column with ones for intercept if included in the model
  if(intrcpt){
    Xc <- cbind(X,rep(1,n))
    unpen<-c(unpen,p+1) #add index of last column for intercept to unpenalised covariates 
  } else{
    Xc <- X
  }
  if(model%in%c("logistic","cox")){
    Xcinit<-Xc
  }
  
  #-3.1.2 Variables used in initialisation of beta (and afterwards) #########################################
  muhat <- array(mu,c(sum(G),nIt+1)); #group means, optional: fixed at mu if given
  tauhatold <- array(tausq,c(sum(G),nIt+1)) #group variances before truncating negative tau to 0, optional: fixed at tausq if given (tausq 1 value for all groups)
  tauhat <- array(max(0,tausq),c(sum(G),nIt+1)) #group variances truncated at 0 for negative values, optional: fixed at tausq if given
  colnames(muhat)<-paste("Itr",0:nIt,sep="")
  colnames(tauhatold)<-paste("Itr",0:nIt,sep="")
  colnames(tauhat)<-paste("Itr",0:nIt,sep="")
  tempRow <- unlist(lapply(1:length(G),function(x){paste("Grouping",x,".G",1:G[x],sep="")}))
  rownames(muhat)<-tempRow;  rownames(tauhatold)<-tempRow;  rownames(tauhat)<-tempRow
  
  weightsMu <- array(NaN,c(sum(G),nIt+1))
  weightsTau <- array(NaN,c(sum(G),nIt+1))
  if(is.nan(mu)){
    partWeightsMu <- array(1,c(m,nIt+1))
    partWeightsMuG<- array(1,c(sum(G),nIt+1))
  }else{
    partWeightsMu <- array(NaN,c(m,nIt+1))
    partWeightsMuG<- array(NaN,c(sum(G),nIt+1))
  }
  if(is.nan(tausq)){
    partWeightsTau <-array(1,c(m,nIt+1))
    partWeightsTauG<-array(1,c(sum(G),nIt+1))
  }else{
    partWeightsTau <-array(NaN,c(m,nIt+1))
    partWeightsTauG<-array(NaN,c(sum(G),nIt+1))
  }
  
  
  #-3.1.3 Variables used in iterations #######################################################################
  if(!missing(X2)){
    if(model=="cox") YpredGR <- array(NaN,c(n2,nIt+1))
    else YpredGR <- array(NaN,c(n2,nIt+1))
    MSEecpc<-rep(NaN,nIt+1)
  } else { 
      YpredGR<-NaN
      MSEecpc<-NaN
      if(!is.nan(compare)){
        Ypredridge<-NaN
        MSEridge<-NaN
      }
  }
  ind0 <- c(); #keep track of index of groups which have 0 variance
  indnot0 <- 1:sum(G) #and keep track of index of groups which have variance > 0
  lambdashat<-array(NaN,c(2,nIt+1,m)) #hyperpenalties for extra level of shrinkage
  colnames(lambdashat)<-paste("Itr",0:nIt,sep="")
  rownames(lambdashat)<-c("PenMu","PenTau")
  lambdashat[1,,]<-hyperlambdas[1]; lambdashat[2,,]<-hyperlambdas[2] #optional: fix extra hyperpenalty if given 
  
  #-3.2 Initial tau and beta ========================================================================================
  if(!silent) print(paste("Estimate global tau^2 (equiv. global ridge penalty lambda)"))
  intrcptGLM <- intrcpt
  MoMinit <- F
  if(grepl("MoM",lambda)){
    lambda<-"CV" #use CV for first betasinit
    intrcptMoM <- intrcpt #memory for next iteration
    MoMinit <- T #use one MoM iteration to update overall tau from initial computed tau
  } 
  #inital tau given
  if(!is.nan(tausq)){
    lambda <- 1/tausq
    if(model=="linear") lambda <- sigmasq/tausq
    if(!is.nan(compare) & compare){ #compare not false
      lambdaridge <- 1/tausq
    }
  }
  if(is.numeric(lambda) & compare) lambdaridge <- lambda

  switch(model,
         'linear'={
           #Use Cross-validation to compute initial lambda (tausq)
           if((!is.nan(compare) & grepl("CV",compare)) | grepl("CV",lambda)){
             #use glmnet to do CV; computationally more expensive but other optimising criteria possible
             if(grepl("glmnet",lambda)){ 
               lambdaGLM<-cv.glmnet(X,Y,nfolds=fold,alpha=0,family=fml,
                                    standardize = F,intercept=intrcpt,
                                    penalty.factor=penfctr,keep=T) #alpha=0 for ridge
             }else{ #use penalized to do CV
               Xsvd <- svd(X[,penfctr!=0])
               XF <- X[,penfctr!=0]%*% Xsvd$v
               if(intrcpt){
                 Xunpen <- cbind(X[,penfctr==0],rep(1,n))
               }else{
                 Xunpen <- cbind(X[,penfctr==0]) #if empty vector, no unpenalised and no intercept
               }
               ol1 <- optL2(Y,penalized=XF, unpenalized =Xunpen,fold=fold,trace=F)
               #ol2 <- optL2(Y,penalized=X[,penfctr!=0],unpenalized=Xunpen, fold=ol1$fold ) #gives same result, but the first is much faster for large p
               itr2<-1
               while(ol1$lambda>10^12 & itr2 < 10){
                 ol1 <- optL2(Y,penalized=XF, unpenalized =Xunpen,fold=fold,trace=F)
                 itr2 <- itr2 + 1
               } 
               if(itr2==10 & ol1$lambda>10^12){
                 ol1$lambda <- 10^12
                 warning("Cross-validated global penalty lambda was >10^12 and set to 10^12")
               }
             }
             
             if((!is.nan(compare) & grepl("CV",compare)) | (!is.nan(compare) & compare==T)){
               # lambdaridge<-lambdaGLM$lambda.min/sd_y*n #using glmnet
               if(grepl("glmnet",lambda)) lambdaridge <- lambdaGLM$lambda.min/sd_y*n #fitted lambda
               else lambdaridge <- ol1$lambda
             } 
             if(grepl("CV",lambda)){
               if(grepl("glmnet",lambda)) lambda <- lambdaGLM$lambda.min/sd_y*n #using glmnet
               else lambda <- ol1$lambda #using penalized
               sigmahat <- sigmasq
               muhat[,1] <- mutrgt
               tauhat[,1] <- 1/lambda
               tautrgt<- 1/lambda
             } 
           }
           #Use ML for sigma estimation and/or initial lambda (tausq) estimate and/or mu
           if(is.nan(sigmasq) | (!is.nan(compare) & grepl("ML",compare)) | grepl("ML",lambda) | is.nan(mutrgt)){
             #Estimate sigma^2, lambda and initial estimate for tau^2 (all betas in one group), mu=0 by default
             #NOTE, TD: not yet possible to include unpenalised covariates in MML
             if(grepl("ML",lambda)){ lambda <- NaN}
             par <- mlestlin(Y,X,lambda=lambda,sigmasq=sigmasq,mu=mutrgt,tausq=tausq) #use maximum marginal likelihood
             lambda <- par[1] 
             sigmahat <- par[2] #sigma could be optimised with CV in the end if not known
             muhat[,1] <- par[3] 
             tauhat[,1] <- par[4]
             
             tautrgt<- par[4] #set target group variance (overall variance if all covariates in one group)
             mutrgt <- par[3] #set target group mean (overall mean if all covariates in one group), 0 by default
             
             if((!is.nan(compare) & grepl("ML",compare)) | (!is.nan(compare) & compare==T)) lambdaridge<- par[1]
           }
           
           #initial estimate for beta
           lambdap <- rep(lambda,p) #px1 vector with penalty for each beta_k, k=1,..,p
           muinitp <- as.vector(c(muhat[,1])%*%Zt) #px1 vector with estimated prior mean for beta_k, k=1,..,p (0 for unpenalised covariates) 
           glmGRtrgt <- glmnet(X,Y,alpha=0,
                               lambda = lambda/n*sd_y,family=fml,
                               offset = X[,!((1:p)%in%unpen)] %*% muinitp, intercept = intrcpt, standardize = F,
                               penalty.factor=penfctr,thresh = 10^-10)
           betasinit <- as.vector(glmGRtrgt$beta)
           betasinit[!((1:p)%in%unpen)] <- betasinit[!((1:p)%in%unpen)] + muinitp
           intrcptinit <- glmGRtrgt$a0
         },
         'logistic'={
           #Use Cross-validation to compute initial lambda (tausq)
           if((!is.nan(compare) & grepl("CV",compare)) | grepl("CV",lambda)){
             #use glmnet to do CV; computationally more expensive but other optimising criteria possible
             if(grepl("glmnet",lambda)){ 
               lambdaGLM<-cv.glmnet(X,Y,nfolds=fold,alpha=0,family=fml,
                                    standardize = F,intercept=intrcpt,
                                    penalty.factor=penfctr,keep=T) #alpha=0 for ridge
             }else{ #use penalized to do CV
               Xsvd <- svd(X[,penfctr!=0])
               XF <- X[,penfctr!=0]%*% Xsvd$v
               if(intrcpt){
                 Xunpen <- cbind(X[,penfctr==0],rep(1,n))
               }else{
                 Xunpen <- cbind(X[,penfctr==0]) #if empty vector, no unpenalised and no intercept
               }
               
               ol1 <- optL2(Y,penalized=XF, unpenalized =Xunpen,fold=fold,trace=F)
               #ol2 <- optL2(Y,penalized=X[,penfctr!=0],unpenalized=Xunpen, fold=ol1$fold ) #gives same result, but the first is much faster for large p
               itr2<-1
               while(ol1$lambda>10^12 & itr2 < 10){
                 ol1 <- optL2(Y,penalized=XF, unpenalized =Xunpen,fold=fold,trace=F)
                 itr2 <- itr2 + 1
                 if(ol1$lambda>10^12){
                   ol1$lambda <- 10^2
                   warning("Cross-validated global penalty lambda was >10^12 and set to 100")
                 }
               } 
               if(itr2==10 & ol1$lambda>10^12){
                 ol1$lambda <- 10^12
                 warning("Cross-validated global penalty lambda was >10^12 and set to 10^12")
               }
               
             }
             if((!is.nan(compare) & grepl("CV",compare)) | (!is.nan(compare) & compare==T)){
               if(grepl("glmnet",lambda)) lambdaridge <- lambdaGLM$lambda.min/sd_y*n #fitted lambda
               else lambdaridge <- ol1$lambda
             } 
             if(grepl("CV",lambda)){
               if(grepl("glmnet",lambda)) lambda <- lambdaGLM$lambda.min/sd_y*n #using glmnet
               else lambda <- ol1$lambda #using penalized
             } 
           }
           #Use approximate ML criterium initial lambda (tausq) estimate 
           if((!is.nan(compare) & grepl("ML",compare)) | grepl("ML",lambda)){
             #Approximate derivative dML/dtausq
             dlMLdtau <- function(param,svdX,Y,UU,rankXXt){
               tausq <- param[1]
               p0 <- param[2]
               
               dlMLdtau <- sum(sapply(1:rankXXt,function(j){
                 T1 <- c(Y-p0)%*%UU[[j]]%*%c(Y-p0) *svdX$d[j]^2
                 T2<- p0*(1-p0)*svdX$d[j]^2*(p0*(1-p0)*svdX$d[j]^2*tausq + 1)
                 T3 <- (p0*(1-p0)*svdX$d[j]^2*tausq + 1)^2
                 return(0.5*(T1-T2)/T3)
               })
               )
               return(dlMLdtau)
             }
             #some variables needed
             svdX<-svd(X)
             svdX$dmin <- 1/svdX$d
             svdX$dmin[svdX$d < 10^-12] <- 0
             rankXXt <- length(svdX$dmin)-sum(svdX$dmin==0)
             U1<-svdX$u[,1:rankXXt]
             UU <- lapply(1:rankXXt,function(j){
               U1[,j]%*%t(U1[,j])
             }
             )
             #fix intercept
             p0hat <- sum(Y)/length(Y)
             
             #find tausq that minimises dlML/dtausq
             opt <- optim((1+p0hat)/(p*p0hat^3),function(x){dlMLdtau(c(x,p0hat),svdX,Y,UU,rankXXt)})
             if( (!is.nan(compare) & grepl("ML",compare)) | (!is.nan(compare) & compare==T)) lambdaridge <- 1/opt$par
             if(grepl("ML",lambda)) lambda <- 1/opt$par
           }
           tauhat[,1] <- 1/lambda
           tautrgt <- 1/lambda
           sigmahat <- 1 #sigma not in model for logistic: set to 1
           muhat[,1] <- mu #use initial mean 0 in logistic setting
           mutrgt <- mutrgt #default: 0
      
           #initial estimate for beta
           lambdap <- rep(lambda,p) #px1 vector with penalty for each beta_k, k=1,..,p
           muinitp<- as.vector(c(muhat[,1])%*%Zt) #px1 vector with estimated prior mean for beta_k, k=1,..,p 
           glmGRtrgt <- glmnet(X,Y,alpha=0,
                               lambda = lambda/n*sd_y,family=fml,
                               offset = X[,!((1:p)%in%unpen)] %*% muinitp, intercept = intrcpt, standardize = F,
                               penalty.factor=penfctr,thresh=10^-20)
           betasinit <- as.vector(glmGRtrgt$beta)
           betasinit[!((1:p)%in%unpen)] <- betasinit[!((1:p)%in%unpen)] + muinitp
           intrcptinit <- glmGRtrgt$a0
         },
         'cox'={
           #Cross-validation lambda
           if((!is.nan(compare) & grepl("CV",compare)) | grepl("CV",lambda)){
             if(grepl("glmnet",lambda)){ 
               lambdaGLM<-cv.glmnet(X,Y,nfolds=fold,alpha=0,family=fml,
                                    standardize = F,
                                    penalty.factor=penfctr,keep=T) #alpha=0 for ridge
             }else{ #use penalized to do CV
               Xsvd <- svd(X[,penfctr!=0])
               XF <- X[,penfctr!=0]%*% Xsvd$v
               Xunpen <- cbind(X[,penfctr==0]) #if empty vector, no unpenalised and no intercept
               
               ol1 <- optL2(Surv(Y[,1],Y[,2]),penalized=XF, unpenalized =Xunpen,fold=fold,trace=F)
               #ol2 <- optL2(Y,penalized=X[,penfctr!=0],unpenalized=Xunpen, fold=ol1$fold ) #gives same result, but the first is much faster for large p
             }
             if((!is.nan(compare) & grepl("CV",compare))| (!is.nan(compare) & compare==T)){
               if(grepl("glmnet",lambda)) lambdaridge <- lambdaGLM$lambda.min/sd_y*n #fitted lambda
               else lambdaridge <- ol1$lambda
             } 
             if(grepl("CV",lambda)){
               if(grepl("glmnet",lambda)) lambda <- lambdaGLM$lambda.min/sd_y*n #fitted lambda
               else lambda <- ol1$lambda
             } 
           }
           sigmahat <- 1 #sigma not in model for cox: set to 1
           muhat[,1] <- 0 #use initial mean 0 in cox setting
           tauhat[,1] <- 1/lambda
           mutrgt <- 0
           tautrgt <- 1/lambda
           
           #initial estimate for beta
           lambdap <- rep(lambda,p) #px1 vector with penalty for each beta_k, k=1,..,p
           muinitp<- as.vector(c(muhat[,1])%*%Zt) #px1 vector with estimated prior mean for beta_k, k=1,..,p
           glmGRtrgt <- glmnet(X,Y,alpha=0,
                               lambda = lambda/n*sd_y,family=fml,
                               offset = X[,!((1:p)%in%unpen)] %*% muinitp, standardize = F,
                               penalty.factor=penfctr,thresh = 10^-10)
           betasinit <- as.vector(glmGRtrgt$beta)
           betasinit[!((1:p)%in%unpen)] <- betasinit[!((1:p)%in%unpen)] + muinitp
           intrcptinit <- glmGRtrgt$a0
         }
  )
  
  #-3.3 Start iterations ========================================================================================
  Itr<-1
  while(Itr<=nIt){
    #-3.3.1 Compute penalty matrix and weight matrix for logistic #############################################
    #copy penalty parameter matrix ridge: add 0 for unpenalised intercept if included
    if(intrcpt | intrcptGLM){
      #Deltac <- diag(c(lambdap,0))
      Deltac <- sparseMatrix(i=1:(length(lambdap)+1),j=1:(length(lambdap)+1),x=c(lambdap,0))
      if(model=="logistic"){
        Deltac<-2*Deltac
        #reweight Xc for logistic model
        expminXb<-exp(-Xcinit%*%c(betasinit,intrcptinit))
        Pinit<-1/(1+expminXb)
        W<-diag(c(sqrt(Pinit*(1-Pinit))))
        Xc<-W%*%Xcinit
      }
    } else{
      #Deltac <- diag(c(lambdap))
      Deltac <- sparseMatrix(i=1:length(lambdap),j=1:length(lambdap),x=c(lambdap))
      if(model=="logistic"){
        Deltac<-2*Deltac
        #reweight Xc for logistic model
        expminXb<-exp(-Xcinit%*%c(betasinit))
        Pinit<-1/(1+expminXb)
        W<-diag(c(sqrt(Pinit*(1-Pinit))))
        Xc<-W%*%Xcinit
      }
      if(model=="cox"){
        Deltac<-2*Deltac
        #reweight Xc for cox model
        expXb<-exp(Xcinit%*%c(betasinit))
        h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXb[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
        H0 <- sapply(Y[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
        
        W <- diag(c(sqrt(H0*expXb)))
        Xc<-W%*%Xcinit
      }
    }
    if(model%in%c("logistic","cox") && all(W==0)){
      #browser()
      if(!silent) print("Overfitting: only 0 in weight matrix W")
      if(!silent) print(paste("Iterating stopped after",Itr-1,"iterations",sep=" "))
      break;
    }
    
    #NOTE: in glmnet not yet unpenalised covariates other than intercept
    
    #-3.3.2 Compute matrices needed for MoM ###################################################################
    # XtXD <- t(Xc)%*%Xc+Deltac
    # XtXDinv <- solve(XtXD)
    # L<-XtXDinv %*% t(Xc)
    pen <- (1:dim(Xc)[2])[!(1:dim(Xc)[2]%in%unpen)] #index covariates to be penalised
    if(p>n){
      if(length(unpen)>0){
        P1<-diag(1,n) - Xc[,unpen]%*%solve(t(Xc[,unpen])%*%Xc[,unpen],t(Xc[,unpen])) #compute orthogonal projection matrix
        eigP1<-eigen(P1,symmetric = T)
        if(!all(round(eigP1$values,digits=0)%in%c(0,1))){
          warning("Check unpenalised covariates")
        }
        eigP1$values<-pmax(0,eigP1$values) #set eigenvalues that are small negative due to numerical errors to 0
        CP1 <- eigP1$vectors %*% diag(sqrt(eigP1$values))
        Xpen <- as.matrix(t(CP1)%*%Xc[,pen]%*%sparseMatrix(i=1:length(pen),j=1:length(pen),x=diag(Deltac)[pen]^(-0.5)))
      } else{
        pen <- 1:p
        CP1<- diag(1,n)
        Xpen <- as.matrix(Xc[,pen]%*%sparseMatrix(i=1:length(pen),j=1:length(pen),x=diag(Deltac)[pen]^(-0.5)))
      }
      
      svdX<-svd(Xpen) #X=UDV^T=RV^T
      svdXR<-svdX$u%*%diag(svdX$d) #R=UD
      L2<-as.matrix(sparseMatrix(i=1:length(pen),j=1:length(pen),x=diag(Deltac)[pen]^(-0.5))%*%svdX$v%*%solve(t(svdXR)%*%svdXR+diag(1,n),t(svdXR)%*%t(CP1)))
      L<-array(0,c(p+intrcpt,n))
      L[pen,]<-L2 #compute only elements corresponding to penalised covariates
      
      R<-Xc
      V2<-sigmahat*apply(L2,1,function(x){sum(x^2)})
      #V2<-apply(L2,1,function(x){sum(x^2)})
      V<-rep(NaN,p+intrcpt)
      V[pen]<-V2 #variance beta ridge estimator
      
      # #should be same as:
      # XtXD <- t(Xc)%*%Xc+Deltac
      # XtXDinv <- solve(XtXD) #inverting pxp matrix really slow, use SVD instead
      # L1<-XtXDinv %*% t(Xc)
      # R<-Xc
      # V1<-sigmahat*apply(L,1,function(x){sum(x^2)})
      # #same as: V <- sigmahat*diag(L%*%R %*% XtXDinv) 
    }else{ #n>p
      XtXD <- as.matrix(t(Xc)%*%Xc+Deltac)
      XtXDinv <- solve(XtXD) #inverting pxp matrix
      L<-XtXDinv %*% t(Xc)
      R<-Xc
      #C<-L%*%R
      V<-sigmahat*apply(L,1,function(x){sum(x^2)})
      #same as: V3 <- sigmahat*diag(L%*%R %*% XtXDinv)
    }
    
    #if lambda=="MoM" inserted in the function, need to update initial beta first before computing group weights
    if(Itr==1 && MoMinit){
      #Compute targets
      if(!is.nan(mu)){mutrgt<-mu}
      if(is.nan(mutrgt)){
        Gamma1<-sum(c(apply(L[pen,],2,sum))*R[,pen])/length(pen)
        Bmu1 <- sum(betasinit[pen]-muinitp[pen]+L[pen,]%*%(R[,pen]%*%muinitp[pen]))/length(pen)
        mutrgt<-Bmu1/Gamma1
      }
      mukhat1<- muinitp[pen] + L[pen,]%*%(R[,pen]%*%(mutrgt-muinitp[pen]))

      #update tau overall
      Btau1 <- sum(pmax((betasinit[pen]^2-mukhat1[pen]^2)/V[pen]-1,0)) / length(pen)
      A1<-sum((t(L[pen,]/c(V[pen]))%*%L[pen,])*(R[,pen]%*%t(R[,pen])))/length(pen)
      tautrgt<-Btau1/A1
      
      lambda <- 1/tautrgt
      tauhat[,1] <- tautrgt
      sigmahat <- 1 #sigma not in model for logistic: set to 1
      muhat[,1] <- mutrgt #use initial mean 0 in logistic setting
      
      if((!is.nan(compare) & grepl("MoM",compare))| (!is.nan(compare) & compare==T)) lambdaridge<-1/tautrgt
      
      #Update beta
      lambdap <- rep(lambda,p) #px1 vector with penalty for each beta_k, k=1,..,p
      muinitp<- as.vector(c(muhat[,1])%*%Zt) #px1 vector with estimated prior mean for beta_k, k=1,..,p 
      if(model=="cox"){
        glmGRtrgt <- glmnet(X,Y,alpha=0,
                        lambda = lambda/n*sd_y,family=fml,
                        offset = X[,!((1:p)%in%unpen)] %*% muinitp, standardize = F,
                        penalty.factor=penfctr)
      }else if(model=="logistic"){
        intrcpt <- intrcptMoM #reset intercept instead of using logit(p0hat)
        glmGRtrgt <- glmnet(X,Y,alpha=0,
                        lambda = lambda/n*sd_y,family=fml,
                        offset = X[,!((1:p)%in%unpen)] %*% muinitp, intercept = intrcpt, standardize = F,
                        penalty.factor=penfctr)
      }else{
        glmGRtrgt <- glmnet(X,Y,alpha=0,
                            lambda = lambda/n/2*sd_y,family=fml,
                            offset = X[,!((1:p)%in%unpen)] %*% muinitp, intercept = intrcpt, standardize = F,
                            penalty.factor=penfctr)
      }
      betasinit <- as.vector(glmGRtrgt$beta)
      betasinit[pen] <- betasinit[pen] + muinitp
      intrcptinit <- glmGRtrgt$a0

      MoMinit <- F
      
      #repeat iteration starting from MoM estimated overall penalty
      next
    }
 
    #-3.3.3 Update group parameters ###########################################################################
    if(nIt>1){
      if(!silent) print(paste("Compute group penalty estimates, iteration ",Itr,"out of maximal ",nIt," iterations."))
    }
    ### Function Method of Moments to compute group weights for (possibly) multiple parameters 
    MoM <- function(Partitions,hypershrinkage=NaN,groupings.grouplvl=NaN,
                    fixWeightsMu=NaN,fixWeightsTau=NaN,pars){
      #Partitions: vector with index of partitions
      #fixWeightsMu,fixWeightsTau: when fixed group weights for different partitions/co-data are given, 
      #                            the MoM-function will calculate partition/co-data weights (without extra shrinkage)
      
      #extract parts of global variables for local copy
      if(length(Partitions)<m | m==1){
        Zt <- Zt[unlist(indGrpsGlobal[Partitions]),,drop=F]
        if(missing(pars)){
          ind0 <- which(unlist(indGrpsGlobal[Partitions])%in%ind0)
          indnot0 <- which(unlist(indGrpsGlobal[Partitions])%in%indnot0)
        }else{
          indnot0 <- pars[[1]]
          ind0 <- pars[[2]]
        }
        if(length(dim(Wminhalf))>1){
          Wminhalf <- Wminhalf[unlist(indGrpsGlobal[Partitions]),unlist(indGrpsGlobal[Partitions]),drop=F]
          # inittauhat <- as.vector(ind[unlist(indGrpsGlobal[Partitions]),]%*%betasinit^2 / 
          #                           apply(ind[unlist(indGrpsGlobal[Partitions]),],1,sum))
          inittauhat <- rep(1,length(indnot0))
        }else inittauhat <- tautrgt  
        G <- G[Partitions] #number of groups in these partitions
        PenGrps <- PenGrps[unlist(indGrpsGlobal[Partitions]),unlist(indGrpsGlobal[Partitions]),drop=F]
        eqfun <- function(gamma,b,A,lam)  return(sum(t(Zt[indnot0,])%*%gamma)/length(pen) ) #equality constraint for average prior variance
      }
      #keep local copies of variables to return
      muhat <- muhat[unlist(indGrpsGlobal[Partitions]),Itr]
      tauhatold <- rep(0,sum(G))
      tauhat <- tauhat[unlist(indGrpsGlobal[Partitions]),Itr]
      lambdashat<-lambdashat[,Itr+1,Partitions] 
      
      #if two shrinkage penalties are given, first is used to select groups, second to shrink estimates
      if(!is.nan(hypershrinkage)){
        temp <- strsplit(hypershrinkage,",") 
        hypershrinkage <- temp[[1]][1]  
        ExtraShrinkage2 <- temp[[1]][-1]
        if(grepl("none",hypershrinkage)){
          if(length(Partitions)==1){
            if(!silent) print(paste("Grouping ",Partitions,": estimate group weights, hypershrinkage type: ",hypershrinkage,sep=""))
          }
        }else{
          if(!silent) print(paste("Grouping ",Partitions,": estimate hyperlambda for ",hypershrinkage," hypershrinkage",sep=""))
        }
      }
      if(length(G)==1 && G==1){
        lambdashat <- c(0,0)
        muhat <- mutrgt
        weightsMu <- NaN
        if(is.nan(tausq)){
          tauhat <- tautrgt
          tauhatold <- tauhat
          weightsTau <- 1
        }
      }else if(!grepl("none",hypershrinkage) & !all(G==1) & length(indnot0)>1){ 
        #-3.3.3|1 With extra shrinkage -----------------------------------------------------------------
        # Use splits to penalise for too many groups
        # Minimise RSS over lambda1, lambda2 to find optimal penalties for shrinkage on group level
        # Splits group randomly in half for nsplits times, INDin: one half of the split
        INDin <- lapply(1:m,function(prt){ 
          if(!(prt%in%Partitions)){return(NULL)}else{
            replicate(nsplits,lapply(groupings[[prt]],function(x){sample(x,floor(length(x)/2),replace=F)}),simplify=F)  
          }
        }) #keep list of m elements such that index same as groupings
        INDout <- lapply(1:m,function(i){ #for each partition
          if(!(i%in%Partitions)){return(NULL)}else{
          lapply(INDin[[i]],function(indin){ #for each split
            lapply(1:length(groupings[[i]]),function(x){groupings[[i]][[x]][!(groupings[[i]][[x]]%in%indin[[x]])]})})
          }
        })
        #INDin[[Partitions]] <- lapply(groupings[Partitions],function(prt){
        #  replicate(nsplits,lapply(prt,function(x){sample(x,floor(length(x)/2),replace=F)}),simplify=F)
        #})
        #INDout <- lapply(Partitions,function(i){ #for each partition
        #  lapply(INDin[[i]],function(indin){ #for each split
        #    lapply(1:G[i],function(x){groupings[[i]][[x]][!(groupings[[i]][[x]]%in%indin[[x]])]})})
        #})
        
        #-3.3.3|1.1 EB estimate group means ============================================================
        muhatp <-as.vector(rep(mu,sum(G))%*%Zt) #px1 vector with estimated prior mean for beta_k, k=1,..,p
        weightsMu <- rep(NaN,sum(G))
        if(is.nan(mu)){
          if(is.nan(lambdashat[1])){
            #-3.3.3|1.1.1 Compute linear system for whole partition ####################################
            Gamma <- matrix(unlist(
              lapply(Partitions,function(i){ #for each partition
                sapply(1:length(Kg[[i]]),function(j){ #for each group
                  #compute row with gamma_{xy}
                  x<-groupings[[i]][[j]]
                  unlist(sapply(Partitions,function(prt){sapply(groupings[[prt]],function(y){sum(L[x,]%*%t(t(R[,y])/Ik[[prt]][y]))/Kg[[i]][j]})}))
                }, simplify="array")
              })
            ),c(sum(G),sum(G)),byrow=T) #reshape to matrix of size sum(G)xsum(G)
            Bmu <- unlist(
              lapply(Partitions,function(i){ #for each partition
                sapply(1:length(Kg[[i]]),function(j){ #for each group
                  x<-groupings[[i]][[j]]
                  sum(betasinit[x]-muinitp[x]+L[x,]%*%(R[,pen]%*%muinitp[pen]))/Kg[[i]][j]
                })
              })
            )
            
            sdGamma <- c(apply(Gamma,2,function(x){sd(x,na.rm=T)}))
            Gamma<-Gamma%*% diag(1/sdGamma) #normalise columns
            
            #-3.3.3|1.1.2 For each split, compute linear system ########################################
            mutrgtG <- mutrgt
            if(length(mutrgt)==1){ mutrgtG <- rep(mutrgt,sum(G))}
            mutrgtG<-diag(sdGamma)%*%mutrgtG
            
            #in-part
            Gammain <- lapply(1:nsplits,function(split){
              matrix(unlist(
                lapply(Partitions,function(i){ #for each partition
                  sapply(1:length(Kg[[i]]),function(j){ #for each group
                    #compute row with gamma_{xy}
                    x<-INDin[[i]][[split]][[j]]
                    #compute row with gamma_{xy}
                    unlist(sapply(Partitions,function(prt){sapply(groupings[[prt]],function(y){sum(L[x,]%*%t(t(R[,y])/Ik[[prt]][y]))/Kg[[i]][j]})}))
                  }, simplify="array")
                })
              ),c(sum(G),sum(G)),byrow=T) %*% diag(1/sdGamma) #reshape to matrix of size sum(G)xsum(G)
            })
            #rhs vector
            Bmuin <- lapply(1:nsplits,function(split){unlist(
              lapply(Partitions,function(i){ #for each partition
                sapply(1:length(Kg[[i]]),function(j){ #for each group
                  x<-INDin[[i]][[split]][[j]]
                  sum(betasinit[x]-muinitp[x]+L[x,]%*%(R[,pen]%*%muinitp[pen]))/Kg[[i]][j]
                })
              })
            )
            })
            
            #weight matrix
            GammainAcc <- lapply(1:nsplits,function(i){
              GammainAcc <- Gammain[[i]] %*% Wminhalf #weight matrix 
            })
            
            #out-part: use Gamma_{out}=Gamma-Gamma_{in}, B_{out}=B-B_{in}
            Gammaout <- lapply(1:nsplits,function(split){
              Gamma-Gammain[[split]]
            })
            Bmuout <- lapply(1:nsplits,function(split){
              Bmu-Bmuin[[split]]
            })
            
            
            #-3.3.3|1.1.3 Define function RSSlambdamu, ################################################
            # using the extra shrinkage penalty function corresponding to parameter hypershrinkage
            rangelambda1 <- c(-100,100)
            switch(hypershrinkage,
                   "ridge"={
                     #standard deviation needed for glmnet
                     sd_Bmuin<- lapply(1:nsplits,function(i){
                       if(length(ind0)>0){
                         sd_Bmuin <- sqrt(var(Bmuin[[i]][indnot0]- 
                                                as.matrix(Gammain[[i]][indnot0,ind0],c(length(indnot0,ind0)))%*%muhat[ind0] -
                                                Gammain[[i]][indnot0,indnot0] %*% mutrgtG[indnot0])*(length(indnot0)-1)/length(indnot0))[1]
                       }else{
                         sd_Bmuin <- sqrt(var(Bmuin[[i]][indnot0]- 
                                                Gammain[[i]][indnot0,indnot0] %*% mutrgtG[indnot0])*(length(indnot0)-1)/length(indnot0))[1]
                       }
                     })
                     RSSlambdamu <- function(lambda1){
                       #Ridge estimates for given lambda
                       lambda1<-exp(lambda1)
                       
                       ### Estimate group means in-part for given lambda1
                       muhatin <- lapply(1:nsplits,function(i){
                         #ridge estimate for group means
                         muhatin <- rep(NaN,sum(G))
                         muhatin[ind0]<-muhat[ind0] #groups with variance 0 keep same prior parameters
                         if(length(ind0)>0){
                           glmMuin <- glmnet(GammainAcc[[i]][indnot0,indnot0],Bmuin[[i]][indnot0]- 
                                               as.matrix(GammainAcc[[i]][indnot0,ind0],c(length(indnot0),length(ind0)))%*%muhat[ind0],
                                             alpha=0,
                                             lambda = 2*lambda1/length(indnot0)*sd_Bmuin[[i]],family="gaussian",
                                             offset = Gammain[[i]][indnot0,indnot0] %*% mutrgtG[indnot0], intercept = F, standardize = F)
                         }else{
                           glmMuin <- glmnet(GammainAcc[[i]][indnot0,indnot0],Bmuin[[i]][indnot0],
                                             alpha=0,
                                             lambda = 2*lambda1/length(indnot0)*sd_Bmuin[[i]],family="gaussian",
                                             offset = Gammain[[i]][indnot0,indnot0] %*% mutrgtG[indnot0], intercept = F, standardize = F)
                         }
                         
                         muhatin[indnot0] <- Wminhalf[indnot0,indnot0] %*% as.vector(glmMuin$beta) + mutrgtG[indnot0]
                         return(muhatin)
                       }
                       ) #group estimate for mu_in
                       
                       ### Compute RSS on left-out part
                       Gammaoutmuin <- lapply(1:nsplits,function(split){Gammaout[[split]][indnot0,]%*%muhatin[[split]]})
                       RSSmu <- sum(sapply(1:nsplits,function(split){sum((Gammaoutmuin[[split]]-Bmuout[[split]][indnot0])^2)/nsplits}))
                       return(RSSmu)
                     }
                   },
                   "lasso"={
                     ### Fit glmnet for global range of lambda
                     fitMu <- lapply(1:nsplits,function(i){
                       if(length(ind0)>0){
                         glmMuin <- glmnet(GammainAcc[[i]][indnot0,indnot0],Bmuin[[i]][indnot0]- 
                                             as.matrix(GammainAcc[[i]][indnot0,ind0],c(length(indnot0),length(ind0)))%*%muhat[ind0],
                                           alpha=1,family="gaussian",
                                           offset = Gammain[[i]][indnot0,indnot0] %*% mutrgtG[indnot0], intercept = F, standardize = F,
                                           thresh = 1e-10)
                       }else{
                         glmMuin <- glmnet(GammainAcc[[i]][indnot0,indnot0],Bmuin[[i]][indnot0],
                                           alpha=1,family="gaussian",
                                           offset = Gammain[[i]][indnot0,indnot0] %*% mutrgtG[indnot0], intercept = F, standardize = F,
                                           thresh = 1e-10)
                       }
                     })
                     
                     RSSlambdamu <- function(lambda1){
                       #Ridge estimates for given lambda
                       lambda1<-exp(lambda1)
                       
                       ### Estimate group means in-part for given lambda1
                       muhatin <- lapply(1:nsplits,function(i){
                         #ridge estimate for group means
                         muhatin <- rep(NaN,sum(G))
                         muhatin[ind0]<-muhat[ind0] #groups with variance 0 keep same prior parameters
                         coefMu<- coef(fitMu[[i]], s = lambda1, exact = F)[-1,]
                         muhatin[indnot0] <- Wminhalf[indnot0,indnot0] %*% as.vector(coefMu) + mutrgtG[indnot0]
                         return(muhatin)
                       }
                       ) #group estimate for mu_in
                       
                       ### Compute RSS on left-out part
                       Gammaoutmuin <- lapply(1:nsplits,function(split){Gammaout[[split]][indnot0,]%*%muhatin[[split]]})
                       RSSmu <- sum(sapply(1:nsplits,function(split){sum((Gammaoutmuin[[split]]-Bmuout[[split]][indnot0])^2)/nsplits}))
                       return(RSSmu)
                     }
                   },
                   "hierLasso"={
                     #TD: acc or not?
                     #Hierarchical overlapping group estimates for given lambda
                     #no target for mu (shrunk to 0)
                     #Gammaxtnd <- lapply(GammainAcc,function(X){return(X[,unlist(groupings.grouplvl)])}) #extend matrix such to create artifical non-overlapping groups
                     Gammaxtnd <- lapply(Gammain,function(X){return(X[,unlist(groupings.grouplvl)])}) #extend matrix such to create artifical non-overlapping groups
                     #create new group indices for Axtnd
                     Kg2 <- c(1,sapply(groupings.grouplvl,length)) #group sizes on group level (1 added to easily compute hier. group numbers)
                     G2 <- length(Kg2)-1
                     groupxtnd <- lapply(2:length(Kg2),function(i){sum(Kg2[1:(i-1)]):(sum(Kg2[1:i])-1)}) #list of indices in each group
                     groupxtnd2 <- unlist(sapply(1:G2,function(x){rep(x,Kg2[x+1])})) #vector with group number
                     
                     ### Fit gglasso for global range of lambda
                     fit1<-lapply(1:nsplits,function(i){
                       gglasso(x=Gammaxtnd[[i]],y=Bmuin[[i]],group = groupxtnd2, loss="ls", 
                               intercept = F, pf = rep(1,G2))
                     })
                     rangelambda1 <- log(range(sapply(fit1,function(i){range(i$lambda)})))
                     
                     RSSlambdamu <- function(lambda1){
                       lambda1<-exp(lambda1)
                       
                       ### Estimate prior taus for given lambda2 (and mutrgt=0)
                       muhatin <- lapply(1:nsplits,function(i){
                         vtilde <- coef(fit1[[i]],s=lambda1)[-1]
                         v<-lapply(groupxtnd,function(g){
                           x<-rep(0,G)
                           x[unlist(groupings.grouplvl)[g]]<-x[unlist(groupings.grouplvl)[g]]+vtilde[g]
                           return(x)
                         })
                         muhatin <- Wminhalf %*% c(apply(array(unlist(v),c(G,G2)),1,sum))
                         return(muhatin)
                       })
                       
                       ### Compute MSE on left-out part
                       Gammaoutmuin <- lapply(1:nsplits,function(split){Gammaout[[split]]%*%muhatin[[split]]})
                       RSSmu <- sum(sapply(1:nsplits,function(split){sum((Gammaoutmuin[[split]]-Bmuout[[split]])^2)/nsplits}))
                       return(RSSmu)
                       
                       # lamb<-seq(exp(rangelambda1[1]),exp(rangelambda1[2]),diff(exp(rangelambda1))/200)
                       # RSS<-sapply(log(lamb),RSSlambdamu)
                       # plot(lamb,RSS)
                     }
                   }
            )
            
            #First find optimal lambda_1
            lambda1<- optimise(RSSlambdamu,rangelambda1)
            lambdashat[1] <- exp(lambda1$minimum)
          }
          
          #-3.3.3|1.1.4 Compute group mean estimates for optimised hyperpenalty lambda #################
          if(lambdashat[1]==0){
            #groups with zero group variance already in muhat
            if(length(ind0)>0){ #only update groups with positive group variance
              muhat[indnot0] <- solve(Gamma[indnot0,indnot0],Bmu[indnot0]-
                                              as.matrix(Gamma[indnot0,ind0],c(length(indnot0),length(ind0)))%*%muhat[ind0])
            }else{
              muhat[indnot0] <- solve(Gamma[indnot0,indnot0],Bmu[indnot0])
            }
            muhat[indnot0] <- diag(1/sdGamma[indnot0]) %*% muhat[indnot0] #restore sd columns mu
          }else{
            #-3.3.3|1.1.5 Compute mu for given hyperpenalty  ###########################################
            switch(hypershrinkage,
                   "ridge"={
                     GammaAcc <- Gamma %*% Wminhalf
                     if(length(ind0)>0){
                       sd_Bmu <- sqrt(var(Bmu[indnot0] - as.matrix(Gamma[indnot0,ind0],c(length(indnot0),length(ind0)))%*%muhat[ind0]
                                          -as.matrix(Gamma[indnot0,indnot0],c(length(indnot0),length(ind0))) %*% mutrgtG[indnot0])*(length(indnot0)-1)/length(indnot0))[1]
                       #ridge estimate for group means
                       glmMu <- glmnet(GammaAcc[indnot0,indnot0],Bmu[indnot0]-
                                         as.matrix(GammaAcc[indnot0,ind0],c(length(indnot0),length(ind0)))%*%muhat[ind0],alpha=0,
                                       lambda = 2*lambdashat[1]/length(indnot0)*sd_Bmu,family="gaussian",
                                       offset = Gamma[indnot0,indnot0] %*% mutrgtG[indnot0], intercept = F, standardize = F)
                     }else{
                       sd_Bmu <- sqrt(var(Bmu[indnot0]
                                          -as.matrix(Gamma[indnot0,indnot0],c(length(indnot0),length(ind0))) %*% mutrgtG[indnot0])*(length(indnot0)-1)/length(indnot0))[1]
                       #ridge estimate for group means
                       glmMu <- glmnet(GammaAcc[indnot0,indnot0],Bmu[indnot0],alpha=0,
                                       lambda = 2*lambdashat[1]/length(indnot0)*sd_Bmu,family="gaussian",
                                       offset = Gamma[indnot0,indnot0] %*% mutrgtG[indnot0], intercept = F, standardize = F)
                     }
                     #groups with variance 0 keep same prior parameters, update other groups
                     muhat[indnot0] <- Wminhalf[indnot0,indnot0] %*% as.vector(glmMu$beta) + mutrgtG[indnot0]
                     muhat[indnot0] <- diag(1/sdGamma[indnot0]) %*% muhat[indnot0] #restore sd columns Gamma
                     
                   },
                   "lasso"={
                     GammaAcc <- Gamma %*% Wminhalf
                     #ridge estimate for group means
                     if(length(ind0)>0){
                       glmMu <- glmnet(GammaAcc[indnot0,indnot0],Bmu[indnot0]-
                                         as.matrix(GammaAcc[indnot0,ind0],c(length(indnot0),length(ind0)))%*%muhat[ind0],
                                       alpha=1,family="gaussian",
                                       offset = Gamma[indnot0,indnot0] %*% mutrgtG[indnot0], intercept = F, standardize = F)
                     }else{
                       glmMu <- glmnet(GammaAcc[indnot0,indnot0],Bmu[indnot0],
                                       alpha=1,family="gaussian",
                                       offset = Gamma[indnot0,indnot0] %*% mutrgtG[indnot0], intercept = F, standardize = F)
                     }
                     coefMu <- coef(glmMu,s=lambdashat[1])
                     #groups with variance 0 keep same prior parameters, update other groups
                     muhat[indnot0] <- Wminhalf[indnot0,indnot0] %*% as.vector(coefMu[-1,]) + mutrgtG[indnot0]
                     muhat[indnot0] <- diag(1/sdGamma[indnot0]) %*% muhat[indnot0] #restore sd columns Gamma
                   },
                   "hierLasso"={
                     #Hierarchical overlapping group estimates for given lambda
                     #no target for mu (shrunk to 0)
                     #GammaAcc <- Gamma %*% Wminhalf
                     #Gammaxtnd <- GammaAcc[,unlist(groupings.grouplvl)] #extend matrix such to create artifical non-overlapping groups
                     Gammaxtnd <- Gamma[,unlist(groupings.grouplvl)] #extend matrix such to create artifical non-overlapping groups
                     #create new group indices for Axtnd
                     Kg2 <- c(1,sapply(groupings.grouplvl,length)) #group sizes on group level (1 added to easily compute hier. group numbers)
                     G2 <- length(Kg2)-1
                     groupxtnd <- lapply(2:length(Kg2),function(i){sum(Kg2[1:(i-1)]):(sum(Kg2[1:i])-1)}) #list of indices in each group
                     groupxtnd2 <- unlist(sapply(1:G2,function(x){rep(x,Kg2[x+1])})) #vector with group number
                     
                     #Hierarchical group lasso estimate for group variances
                     fit1<-gglasso(x=Gammaxtnd,y=Bmu,group = groupxtnd2, loss="ls", 
                                   intercept = F, pf = rep(1,G2),lambda=lambdashat[1])
                     vtilde <- coef(fit1,s=lambdashat[1])[-1]
                     v<-lapply(groupxtnd,function(g){
                       x<-rep(0,G)
                       x[unlist(groupings.grouplvl)[g]]<-x[unlist(groupings.grouplvl)[g]]+vtilde[g]
                       return(x)
                     })
                     muhat <- Wminhalf %*% c(apply(array(unlist(v),c(G,G2)),1,sum))
                     muhat <- diag(1/sdGamma) %*% muhat #restore sd columns A
                   })
          }
          
          weightsMu <- muhat*p/sum(as.vector(c(muhat)%*%Zt))
          muhatp <-as.vector(c(muhat)%*%Zt) #px1 vector with estimated prior mean for beta_k, k=1,..,p
          # if(normalise){ #T by default
          #   C<-mutrgt*p/sum(muhatp)
          #   muhat[,Itr+1]<-muhat[,Itr+1]*C
          #   muhatp <-as.vector(c(muhat[,Itr+1])%*%Zt) #px1 vector with estimated prior mean for beta_k, k=1,..,p
          # }
          
          #should be same as:
          #muhat2 <- solve(t(Gamma)%*%Gamma+lambdashat[1]*diag(weights),t(Gamma)%*%Bmu+lambdashat[1]*diag(weights)%*%rep(mutrgt,G))
          #muhat2 <- solve(t(Gamma)%*%diag(c(Kg))%*%Gamma+lambdashat[1]*diag(1,G),t(Gamma)%*%diag(c(Kg))%*%Bmu+lambdashat[1]*diag(1,G)%*%rep(mutrgt,G))
        }
        
        #-3.3.3|1.2 EB estimate group variances =========================================================
        weightsTau <- rep(1,sum(G))
        if(is.nan(tausq)){
          if(is.nan(lambdashat[2])){
            #-3.3.3|1.2.1 Compute linear system for whole partition #####################################
            Btau <- unlist(
              lapply(Partitions,function(i){ #for each partition
                sapply(1:length(Kg[[i]]),function(j){ #for each group
                  if(j%in%ind0) return(NaN)
                  #compute row with gamma_{xy}
                  x<-groupings[[i]][[j]]
                  sum(pmax(0,(betasinit[x]^2-(muinitp[x]+L[x,]%*%(R[,pen]%*%(muhatp[pen]-muinitp[pen])))^2)/V[x]-1))/Kg[[i]][j]
                  #sum((betasinit[x]^2-(muinitp[x]+L[x,]%*%(R[,pen]%*%(muhatp[pen]-muinitp[pen])))^2)/V[x]-1)/Kg[[i]][j]
                })
              })
            )
            A <- matrix(unlist(
              lapply(Partitions,function(i){ #for each partition
                sapply(1:length(Kg[[i]]),function(j){ #for each group
                  if(j%in%ind0) return(rep(NaN,sum(G)))
                  #compute row with gamma_{xy}
                  x<-groupings[[i]][[j]]
                  #compute row with gamma_{xy}
                  unlist(sapply(Partitions,function(prt){sapply(groupings[[prt]],function(y){
                    sum(t(c(1/V[x])*L[x,])%*%L[x,]*(R[,y]%*%(t(R[,y])/c(Ik[[prt]][y]))))/Kg[[i]][j]
                  })}))
                }, simplify="array")
              })
            ),c(sum(G),sum(G)),byrow=T)*tautrgt #reshape to matrix of size sum(G)xsum(G)
            
            
            #if(Itr==2) browser()
            #-3.3.3|1.2.2 For each split, compute linear system #########################################
            tautrgtG <- tautrgt
            if(length(tautrgt)==1){tautrgtG<-rep(tautrgt,sum(G))}
            tautrgtG[ind0]<-0 
            
            #in-part
            flag <- T; itr2 <- 1
            indNewSplits <- 1:nsplits; 
            Btauin <- list(); Btauout <- list()
            while(flag & itr2 <= 50){
              Btauin[indNewSplits] <- lapply(indNewSplits,function(split){unlist(
                lapply(Partitions,function(i){ #for each partition
                  sapply(1:length(Kg[[i]]),function(j){ #for each group
                    if(j%in%ind0) return(NaN)
                    #compute row with gamma_{xy}
                    x<-INDin[[i]][[split]][[j]]
                    sum(pmax(0,(betasinit[x]^2-(muinitp[x]+L[x,]%*%(R[,pen]%*%(muhatp[pen]-muinitp[pen])))^2)/V[x]-1))/Kg[[i]][j]
                    #sum((betasinit[x]^2-(muinitp[x]+L[x,]%*%(R[,pen]%*%(muhatp[pen]-muinitp[pen])))^2)/V[x]-1)/Kg[[i]][j]
                  })
                })
              )
              })
              
              #check split: at least two elements of Btauin (of selected groups) should be larger than 0
              checkSplit <- sapply(Btauin,function(b){sum(b[!is.nan(b)]!=0)>=2 }) 
              if(all(checkSplit)){ #all splits are fine
                flag <- F
              }else{ 
                itr2 <- itr2 + 1
                indNewSplits <- which(!checkSplit) #index of splits that have to be resampled
                #resample split
                INDin[[Partitions]][indNewSplits] <- replicate(length(indNewSplits),lapply(groupings[[Partitions]],
                                                               function(x){sample(x,floor(length(x)/2),replace=F)}),simplify=F)  
                INDout[[Partitions]][indNewSplits] <- lapply(INDin[[Partitions]][indNewSplits],function(indin){ #for each split
                  lapply(1:length(groupings[[Partitions]]),
                         function(x){groupings[[Partitions]][[x]][!(groupings[[Partitions]][[x]]%in%indin[[x]])]})})
              }
            }
            if(itr2==51) warning("Check splits")
            
            Ain <- lapply(1:nsplits,function(split){
              matrix(unlist(
                lapply(Partitions,function(i){ #for each partition
                  sapply(1:length(Kg[[i]]),function(j){ #for each group
                    if(j%in%ind0) return(rep(NaN,sum(G)))
                    #compute row with gamma_{xy}
                    x<-INDin[[i]][[split]][[j]]
                    #compute row with gamma_{xy}
                    unlist(sapply(Partitions,function(prt){sapply(groupings[[prt]],function(y){
                      sum(t(c(1/V[x])*L[x,])%*%L[x,]*(R[,y]%*%(t(R[,y])/c(Ik[[prt]][y]))))/Kg[[i]][j]
                    })}))
                  }, simplify="array")
                })
              ),c(sum(G),sum(G)),byrow=T) *tautrgt #reshape to matrix of size sum(G)xsum(G)
            })
            #weight matrix
            AinAcc <- lapply(1:nsplits,function(i){
              AinAcc <- Ain[[i]] %*% Wminhalf #weight matrix 
            })
            
            Btauout <- lapply(1:nsplits,function(split){unlist(
              lapply(Partitions,function(i){ #for each partition
                sapply(1:length(Kg[[i]]),function(j){ #for each group
                  if(j%in%ind0) return(NaN)
                  #compute row with gamma_{xy}
                  x<-INDout[[i]][[split]][[j]]
                  sum(pmax(0,(betasinit[x]^2-(muinitp[x]+L[x,]%*%(R[,pen]%*%(muhatp[pen]-muinitp[pen])))^2)/V[x]-1))/Kg[[i]][j]
                  #sum((betasinit[x]^2-(muinitp[x]+L[x,]%*%(R[,pen]%*%(muhatp[pen]-muinitp[pen])))^2)/V[x]-1)/Kg[[i]][j]
                  #sum(pmax(0,(betasinit[x]^2-(muinitp[x]+L[x,]%*%(R[,pen]%*%(muhatp[pen]-muinitp[pen])))^2)-V[x]))/Kg[[i]][j]
                })
              })
            )
            })
            # Btauout <- lapply(1:nsplits,function(split){
            #   Btau - Btauin[[split]]
            # })
            Aout <- lapply(1:nsplits,function(split){
              A - Ain[[split]]
            })
            
            #-3.3.3|1.2.3 Define function RSSlambdatau, #################################################
            # using the extra shrinkage penalty function corresponding to parameter hypershrinkage
            rangelambda2 <- c(10^-9,10^9)
            switch(hypershrinkage,
                   "ridge"={
                     tautrgtG[indnot0] <- 1
                     meanWhalf <- mean(diag(Wminhalf)^-1)
                     #standard deviation needed for glmnet
                     sd_Btauin<- lapply(1:nsplits,function(i){
                       sd_Btauin <- sqrt(var(Btauin[[i]][indnot0] - Ain[[i]][indnot0,indnot0] %*% tautrgtG[indnot0])*(length(indnot0)-1)/length(indnot0))[1]
                     })
                     
                     #function to compute tau for linear system given a hyperpenalty lambda2
                     taus <- function(lambda2){
                       sd_Btau <- sqrt(var(Btau[indnot0] - A[indnot0,indnot0] %*% tautrgtG[indnot0])*(length(indnot0)-1)/length(indnot0))[1]
                       
                       taus <- rep(0,G)
                       Aacc <- A %*% (Wminhalf * meanWhalf)
                       #ridge estimate for group variances
                       glmTau <- glmnet(Aacc[indnot0,indnot0],Btau[indnot0],alpha=0,
                                          lambda = 2*lambda2/length(indnot0)*sd_Btau,family="gaussian",
                                          offset = Aacc[indnot0,indnot0] %*% tautrgtG[indnot0], intercept = F, standardize = F)
                       taus[indnot0] <- (Wminhalf[indnot0,indnot0]*meanWhalf) %*% as.vector(glmTau$beta) + tautrgtG[indnot0] 
                       #taus <- pmax(taus,0) #truncate at 0
                       return(taus)
                     }
                     
                     #function to compute the Residual Sum of Squares on the splits given lambda2
                     RSSlambdatau <- function(lambda2){
                       #Ridge estimates for given lambda
                       ### Estimate prior taus for given lambda2 and mutrgt
                       tauhatin <- lapply(1:nsplits,function(i){
                         tauhatin <- rep(NaN,sum(G))
                         tauhatin[ind0] <- 0
                         #ridge estimate for group variances
                         glmTauin <- glmnet(AinAcc[[i]][indnot0,indnot0],Btauin[[i]][indnot0],alpha=0,
                                            lambda = 2*lambda2/length(indnot0)*sd_Btauin[[i]],family="gaussian",
                                            offset = Ain[[i]][indnot0,indnot0] %*% tautrgtG[indnot0], intercept = F, standardize = F)
                         tauhatin[indnot0] <- Wminhalf[indnot0,indnot0] %*% as.vector(glmTauin$beta) + tautrgtG[indnot0] 
                         #tauhatin <- pmax(tauhatin,0)
                         return(tauhatin)
                       })
                       
                       ### Compute MSE on left-out part
                       Aouttauin <- lapply(1:nsplits,function(split){Aout[[split]][indnot0,]%*%tauhatin[[split]]})
                       RSStau <- sum(sapply(1:nsplits,function(split){sum((Aouttauin[[split]]-Btauout[[split]][indnot0])^2)/nsplits}))
                       return(RSStau)
                     }
                   },
                   "ridgeGAM"={
                     tautrgtG[indnot0] <- 1
                     
                     #function to compute tau for linear system given a hyperpenalty lambda2
                     taus <- function(lambda2){
                       taus <- rep(0,G)
                       
                       #ridge estimate for group variances with fused penalty for overlapping groups
                       dat<-list(Y=Btau[indnot0],X=A[indnot0,indnot0])
                       gamTau <- gam(Y~ 0 + X, data=dat,
                                       family="gaussian",offset = A[indnot0,indnot0] %*% tautrgtG[indnot0],
                                       paraPen=list(X=list(S1=PenGrps,sp=lambda2)))
                       taus[indnot0] <- gamTau$coefficients + tautrgtG[indnot0]
                       
                       return(taus)
                     }
                     
                     #function to compute the Residual Sum of Squares on the splits given lambda2
                     RSSlambdatau <- function(lambda2){

                       ### Estimate prior taus for given lambda2 and mutrgt
                       tauhatin <- lapply(1:nsplits,function(i){
                         tauhatin <- rep(NaN,sum(G))
                         tauhatin[ind0] <- 0
                         #ridge estimate for group variances
                         dat<-list(Y=Btauin[[i]][indnot0],X=Ain[[i]][indnot0,indnot0])
                         gamTauin <- gam(Y~ 0 + X, data=dat,
                                         family="gaussian",offset = Ain[[i]][indnot0,indnot0] %*% tautrgtG[indnot0],
                                         paraPen=list(X=list(S1=PenGrps,sp=lambda2)))
                         tauhatin[indnot0] <- gamTauin$coefficients + tautrgtG[indnot0]
                         
                         return(tauhatin)
                       })
                       
                       ### Compute MSE on left-out part
                       Aouttauin <- lapply(1:nsplits,function(split){Aout[[split]][indnot0,]%*%tauhatin[[split]]})
                       RSStau <- sum(sapply(1:nsplits,function(split){sum((Aouttauin[[split]]-Btauout[[split]][indnot0])^2)/nsplits}))
                       return(RSStau)
                     }
                   },
                   "lasso"={#TD: adapt lasso&hierLasso, check other hyperpenalties on inclusion function taus
                     meanWhalf <- mean(diag(Wminhalf)^-1)
                     ### Fit glmnet lasso for global range of lambda
                     fitTau <- lapply(1:nsplits,function(i){
                       glmTauin <- glmnet(AinAcc[[i]][indnot0,indnot0]*meanWhalf,Btauin[[i]][indnot0],
                                          alpha=1,family="gaussian",
                                          intercept = F, standardize = F,
                                          thresh=1e-6)
                     })
                     
                     #function to compute tau for linear system given a hyperpenalty lambda2
                     taus <- function(lambda2){
                       taus <- rep(0,G)
                       Aacc <- A %*% (Wminhalf * meanWhalf)
                       glmTau <- glmnet(Aacc[indnot0,indnot0],Btau[indnot0],
                                        alpha=1,family="gaussian",
                                        intercept = F, standardize = F)
                       coefTau <- coef(glmTau,s=lambda2,exact=T,
                                       x=Aacc[indnot0,indnot0],y=Btau[indnot0])
                       taus[indnot0] <- (Wminhalf[indnot0,indnot0]*meanWhalf) %*% as.vector(coefTau[-1,])
                       return(taus)
                     }
                     
                     #function to compute the Residual Sum of Squares on the splits given lambda2
                     RSSlambdatau <- function(lambda2){
                       #Ridge estimates for given lambda
                       ### Estimate prior taus for given lambda2 and mutrgt
                       tauhatin <- lapply(1:nsplits,function(i){
                         tauhatin <- rep(NaN,sum(G))
                         tauhatin[ind0] <- 0
                         #ridge estimate for group variances
                         coefTau<- coef(fitTau[[i]], s = lambda2, exact = T,
                                        x=AinAcc[[i]][indnot0,indnot0]*meanWhalf,y=Btauin[[i]][indnot0])[-1,]
                         tauhatin[indnot0] <- (Wminhalf[indnot0,indnot0]*meanWhalf) %*% as.vector(coefTau)
                         return(tauhatin)
                       })
                       
                       ### Compute MSE on left-out part
                       Aouttauin <- lapply(1:nsplits,function(split){Aout[[split]][indnot0,]%*%tauhatin[[split]]})
                       RSStau <- sum(sapply(1:nsplits,function(split){sum((Aouttauin[[split]]-Btauout[[split]][indnot0])^2)/nsplits}))
                       return(RSStau)
                     }
                   },
                   "hierLasso"={
                     #Hierarchical overlapping group estimates for given lambda
                     #no target for tau (shrunk to 0)
                     #remove groups that are already set to 0
                     if(length(groupings.grouplvl)!=length(indnot0)){
                       INDgrps2 <- lapply(groupings.grouplvl[indnot0],function(x){x[x%in%indnot0]})
                     }else{
                       INDgrps2 <- groupings.grouplvl
                     }
                     
                     #Axtnd <- lapply(AinAcc,function(A){return(A[indnot0,unlist(INDgrps2),drop=F])}) #extend matrix such to create artifical non-overlapping groups
                     Axtnd <- lapply(Ain,function(A){return(A[indnot0,unlist(INDgrps2),drop=F])}) #extend matrix such to create artifical non-overlapping groups
                     #create new group indices for Axtnd
                     Kg2 <- c(1,sapply(INDgrps2,length)) #group sizes on group level (1 added to easily compute hier. group numbers)
                     G2 <- length(Kg2)-1
                     groupxtnd <- lapply(2:length(Kg2),function(i){sum(Kg2[1:(i-1)]):(sum(Kg2[1:i])-1)}) #list of indices in each group
                     groupxtnd2 <- unlist(sapply(1:G2,function(x){rep(x,Kg2[x+1])})) #vector with group number
                     
                     ### Fit gglasso for global range of lambda
                     fit2<-lapply(1:nsplits,function(i){
                       gglasso(x=Axtnd[[i]],y=Btauin[[i]][indnot0],group = groupxtnd2, loss="ls",
                               intercept = F, pf = rep(1,G2))
                     })
                     #rangelambda2 <- range(sapply(fit2,function(i){range(i$lambda)}))
                     
                     #Find grid to search optimal lambda over
                     taus <- function(lambda2){
                       taus <- rep(0,G)
                       Axtnd <- A[indnot0,unlist(INDgrps2),drop=F] #extend matrix such to create artifical non-overlapping groups
                       
                       #create new group indices for Axtnd
                       Kg2 <- c(1,sapply(INDgrps2,length)) #group sizes on group level (1 added to easily compute hier. group numbers)
                       G2 <- length(Kg2)-1
                       groupxtnd <- lapply(2:length(Kg2),function(i){sum(Kg2[1:(i-1)]):(sum(Kg2[1:i])-1)}) #list of indices in each group
                       groupxtnd2 <- unlist(sapply(1:G2,function(x){rep(x,Kg2[x+1])})) #vector with group number
                       
                       #Hierarchical group lasso estimate for group variances
                       fit2<-gglasso(x=Axtnd,y=Btau[indnot0],group = groupxtnd2, loss="ls",
                                     intercept = F, pf = rep(1,G2),lambda=lambda2)
                       tauhat <- rep(0,sum(G))
                       vtilde <- coef(fit2,s=lambdashat[2])[-1]
                       v<-lapply(groupxtnd,function(g){
                         x<-rep(0,sum(G))
                         x[unlist(INDgrps2)[g]]<-x[unlist(INDgrps2)[g]]+vtilde[g]
                         return(x)
                       })
                       #tauhatold[indnot0] <- Wminhalf[indnot0,indnot0] %*% c(apply(array(unlist(v),c(sum(G),G2)),1,sum))[indnot0]
                       taus[indnot0] <- c(apply(array(unlist(v),c(sum(G),G2)),1,sum))[indnot0]
                       return(taus)
                     }
                     
                     #function to compute the Residual Sum of Squares on the splits given lambda2
                     RSSlambdatau <- function(lambda2){
                       ### Estimate prior taus for given lambda2 (and mutrgt=0)
                       tauhatin <- lapply(1:nsplits,function(i){
                         tauhatin <- rep(0,sum(G))
                         vtilde <- coef(fit2[[i]],s=lambda2)[-1]
                         v<-lapply(groupxtnd,function(g){
                           x<-rep(0,sum(G))
                           x[unlist(INDgrps2)[g]]<-x[unlist(INDgrps2)[g]]+vtilde[g]
                           return(x)
                         })
                         #tauhatin[indnot0] <- Wminhalf[indnot0,indnot0] %*% c(apply(array(unlist(v),c(sum(G),G2)),1,sum))[indnot0]
                         tauhatin[indnot0] <- c(apply(array(unlist(v),c(sum(G),G2)),1,sum))[indnot0]
                         tauhatin[tauhatin<0] <- 0
                         return(tauhatin)
                       })
                       
                       ### Compute MSE on left-out part
                       Aouttauin <- lapply(1:nsplits,function(split){Aout[[split]][indnot0,indnot0]%*%tauhatin[[split]][indnot0]})
                       RSStau <- sum(sapply(1:nsplits,function(split){sum((Aouttauin[[split]]-Btauout[[split]][indnot0])^2)/nsplits}))
                       return(RSStau)
                     }
                   },
                   "ridge+positive"={
                     meanWhalf <- mean(diag(Wminhalf)^-1)
                     trgt <- 1/diag(Wminhalf)/meanWhalf
                     inittauhat <- diag(Wminhalf)^(-1)/meanWhalf
                     
                     #define function for MSE penalised with ridge prior with target 1
                     penMSE <- function(gamma,b,A,lam){ 
                       return(sum((b-A%*%gamma)^2) + lam*sum((gamma-trgt[indnot0])^2)) }
                     
                     #function to compute tau for linear system given a hyperpenalty lambda2
                     taus <- function(lambda2){
                       taus <- rep(0,G)
                       Aacc <- A %*% (Wminhalf * meanWhalf)
                       initSelected <- inittauhat[indnot0]
                       fitTau <- solnp(par = initSelected, fun=penMSE, b=Btau[indnot0],
                                       A=Aacc[indnot0,indnot0], lam=lambda2,
                                       LB = rep(0,length(indnot0)),control=list(trace=0))
                       taus[indnot0] <- (Wminhalf[indnot0,indnot0] *meanWhalf) %*%as.vector(fitTau$pars)
                       return(taus)
                     }
                     
                     #function to compute the Residual Sum of Squares on the splits given lambda2
                     RSSlambdatau <- function(lambda2){
                       #Ridge estimates for given lambda
                       ### Estimate prior taus for given lambda2 and mutrgt
                       tauhatin <- lapply(1:nsplits,function(i){
                         tauhatin <- rep(NaN,sum(G))
                         tauhatin[ind0] <- 0
             
                         Ainacc <- Ain[[i]]%*%(Wminhalf * meanWhalf)
                         fitTauin <- solnp(par = inittauhat[indnot0], fun=penMSE, b=Btauin[[i]][indnot0],
                                           A=Ainacc[indnot0,indnot0],lam=lambda2,
                                           LB = rep(0,length(indnot0)),control=list(trace=0))
                         tauhatin[indnot0] <- (Wminhalf[indnot0,indnot0] * meanWhalf) %*% as.vector(fitTauin$pars)
                         return(tauhatin)
                       })
                       ### Compute MSE on left-out part
                       Aouttauin <- lapply(1:nsplits,function(split){Aout[[split]][indnot0,]%*%tauhatin[[split]]})
                       RSStau <- sum(sapply(1:nsplits,function(split){sum((Aouttauin[[split]]-Btauout[[split]][indnot0])^2)/nsplits}))
                       return(RSStau)
                     }
                   },
                   "ridgeGAM+positive"={
                     trgt <- 1
                     
                     #define function for MSE penalised with ridge prior with target 1
                     penMSE <- function(gamma,b,A,lam){ 
                       return(sum((b-A%*%gamma)^2) + 
                                lam*sum((gamma-trgt)%*%PenGrps%*%(gamma-trgt))) }
                     
                     taus <- function(lambda2){
                       taus <- rep(0,sum(G))
                       initSelected <- inittauhat[indnot0]
                       fitTau <- solnp(par = initSelected, fun=penMSE, b=Btau[indnot0], 
                                         A=A[indnot0,indnot0], lam=lambda2,
                                         LB = rep(0,length(indnot)),control=list(trace=0))
                       tauhatin[indnot0] <- as.vector(fitTauin$pars)
                       return(tauhatin)
                     }
                     
                     RSSlambdatau <- function(lambda2){
                       #Ridge estimates for given lambda
                       ### Estimate prior taus for given lambda2 and mutrgt
                       tauhatin <- lapply(1:nsplits,function(i){
                         tauhatin <- rep(0,sum(G))
                         
                         initSelected <- inittauhat[indnot0]
                         fitTauin <- solnp(par = initSelected, fun=penMSE, b=Btauin[[i]][indnot0], 
                                           A=Ain[[i]][indnot0,indnot0], lam=lambda2,
                                           LB = rep(0,length(indnot0)),control=list(trace=0))
                         tauhatin[indnot0] <- as.vector(fitTauin$pars)
                         return(tauhatin)
                       })
                       
                       ### Compute MSE on left-out part
                       Aouttauin <- lapply(1:nsplits,function(split){Aout[[split]][indnot0,]%*%tauhatin[[split]]})
                       RSStau <- sum(sapply(1:nsplits,function(split){sum((Aouttauin[[split]]-Btauout[[split]][indnot0])^2)/nsplits}))
                       return(RSStau)
                     }
                   },
                   "invgamma+mean1"={
                     meanWhalf <- mean(diag(Wminhalf)^-1)
                     
                     #define function for MSE penalised with inverse gamma prior with mean 1
                     inittauhat <- diag(Wminhalf)^(-1)/meanWhalf
                     #MSE penalised by inverse gamma penalty
                     penMSE <- function(gamma,b,A,lam){ 
                       Kg <- diag(Wminhalf)^(-2) #group sizes
                       alphaIG <- pmax(1,2 + (lam-1/min(Kg))*Kg) #alpha in range [1,infty)
                       betaIG <- pmax(0,sqrt(Kg)* (1+(lam-1/min(Kg))* Kg) / meanWhalf) #beta in range [0,infty)
                       
                       minlogLikeInvGamma <- (alphaIG[indnot0] + 1)*log(gamma) + betaIG[indnot0]/gamma
                       return(sum((b-A%*%gamma)^2) + sum(minlogLikeInvGamma) ) }
                     
                     #function to compute tau for linear system given a hyperpenalty lambda2
                     taus <- function(lambda2){
                       taus <- rep(0,G)
                       Aacc <- A %*% (Wminhalf * meanWhalf)
                       initSelected <- inittauhat[indnot0]
                       fitTau <- solnp(par = initSelected, fun=penMSE, b=Btau[indnot0],
                                       A=Aacc[indnot0,indnot0], lam=lambda2,
                                       LB = rep(0,length(indnot0)),control=list(trace=0))
                       taus[indnot0] <- (Wminhalf[indnot0,indnot0] *meanWhalf) %*%as.vector(fitTau$pars)
                       return(taus)
                     }
                     
                     #function to compute the Residual Sum of Squares on the splits given lambda2
                     RSSlambdatau <- function(lambda2){
                       ### Estimate prior taus for given lambda2
                       tauhatin <- lapply(1:nsplits,function(i){
                         Ainacc <- Ain[[i]]%*%(Wminhalf * meanWhalf)
                         tauhatin <- rep(NaN,sum(G))
                         tauhatin[ind0] <- 0
                         initSelected <- inittauhat[indnot0]
                         fitTauin <- solnp(par = initSelected, fun=penMSE, b=Btauin[[i]][indnot0],
                                           A=Ainacc[indnot0,indnot0],lam=lambda2,
                                           LB = rep(0,length(indnot0)),control=list(trace=0))
                         tauhatin[indnot0] <- (Wminhalf[indnot0,indnot0] *meanWhalf)%*%as.vector(fitTauin$pars)
                         return(tauhatin)
                       })
                       
                       ### Compute MSE on left-out part
                       Aouttauin <- lapply(1:nsplits,function(split){Aout[[split]][indnot0,]%*%tauhatin[[split]]})
                       RSStau <- sum(sapply(1:nsplits,function(split){sum((Aouttauin[[split]]-Btauout[[split]][indnot0])^2)/nsplits}))
                       return(RSStau)
                     }
                   },
                   "invgamma+mode1"={
                     meanWhalf <- mean(diag(Wminhalf)^-1)
                     
                     #define function for MSE penalised with inverse gamma prior with mean 1
                     inittauhat <- diag(Wminhalf)^(-1)/meanWhalf

                     #MSE penalised by inverse gamma penalty
                     penMSE <- function(gamma,b,A,lam){ 
                       Kg <- diag(Wminhalf)^(-2) #group sizes
                       prmsIG<-prmsIGMode1(lam,Kg)
                       alphaIG <- prmsIG[[1]]
                       betaIG<-prmsIG[[2]] * diag(Wminhalf)^(-1)/meanWhalf
                       
                       minlogLikeInvGamma <- (alphaIG[indnot0] + 1)*log(gamma) + betaIG[indnot0]/gamma
                       return(sum((b-A%*%gamma)^2) + sum(minlogLikeInvGamma) ) }
                    
                     #function to compute tau for linear system given a hyperpenalty lambda2
                     taus <- function(lambda2){
                       taus <- rep(0,G)
                       Aacc <- A %*% (Wminhalf * meanWhalf)
                       initSelected <- inittauhat[indnot0]
                       fitTau <- solnp(par = initSelected, fun=penMSE, b=Btau[indnot0],
                                       A=Aacc[indnot0,indnot0], lam=lambda2,
                                       LB = rep(0,length(indnot0)),control=list(trace=0))
                       taus[indnot0] <- (Wminhalf[indnot0,indnot0] *meanWhalf) %*%as.vector(fitTau$pars)
                       return(taus)
                     }
                    
                     #function to compute the Residual Sum of Squares on the splits given lambda2
                     RSSlambdatau <- function(lambda2){
                       ### Estimate prior taus for given lambda2
                       tauhatin <- lapply(1:nsplits,function(i){
                         Ainacc <- Ain[[i]]%*%(Wminhalf * meanWhalf)
                         tauhatin <- rep(NaN,sum(G))
                         tauhatin[ind0] <- 0
                         initSelected <- inittauhat[indnot0]
                         fitTauin <- solnp(par = initSelected, fun=penMSE, b=Btauin[[i]][indnot0],
                                           A=Ainacc[indnot0,indnot0],lam=lambda2,
                                           LB = rep(0,length(indnot0)),control=list(trace=0))
                         tauhatin[indnot0] <- (Wminhalf[indnot0,indnot0] *meanWhalf)%*%as.vector(fitTauin$pars)
                         return(tauhatin)
                       })
                       
                       ### Compute MSE on left-out part
                       Aouttauin <- lapply(1:nsplits,function(split){Aout[[split]][indnot0,]%*%tauhatin[[split]]})
                       RSStau <- sum(sapply(1:nsplits,function(split){sum((Aouttauin[[split]]-Btauout[[split]][indnot0])^2)/nsplits}))
                       return(RSStau)
                     }
                   },
                   "gamma+positive"={
                     #define function for MSE penalised with gamma prior with mean 1
                     spike<-0.001
                     penMSE <- function(gamma2,b,A,lam){ 
                       #logLikeInvGamma <- sapply(gamma,function(x){logdinvgamma(x,alp=lam,bet=lam-1)})
                       logLikeGamma <- sapply(gamma2,function(x){
                         #return(lam*log(lam)-log(gamma(lam))+(lam-1)*log(x)-lam*x)
                         if(x==0) return(log(spike))
                         return(log(1-spike)+(lam-1)*log(x)-lam*x)
                       })
                       return(sum((b-A%*%gamma2)^2) - sum(logLikeGamma) ) }
                     #eqfun2 <- function(gamma,b,A,lam)  return(sum(t(Zt[indnot0,])%*%(Wminhalf[indnot0,indnot0]%*%gamma))/length(pen) ) #equality constraint for average prior variance
                     
                     RSSlambdatau <- function(lambda2){
                       if(lambda2<100) lambda2<-log(exp(lambda2)+1)

                       browser()
                       ### Estimate prior taus for given lambda2
                       #compute tauhat for first split
                       i<-1
                       #Ainacc <- Ain[[i]]%*%Wminhalf
                       tauhatin1 <- rep(NaN,sum(G))
                       tauhatin1[ind0] <- 0
                       fitTauin1 <- solnp(par = inittauhat, fun=penMSE, b=Btauin[[i]][indnot0], 
                                          A=Ain[[i]][indnot0,indnot0],lam=lambda2,
                                          LB = rep(0,G), eqfun=eqfun, eqB = 1,control=list(trace=0))
                       #tauhatin1[indnot0] <- Wminhalf%*%as.vector(fitTauin1$pars)
                       tauhatin1[indnot0] <- as.vector(fitTauin1$pars)
                       
                       tauhatin <- lapply(2:nsplits,function(i){
                         #Ainacc <- Ain[[i]]%*%Wminhalf
                         tauhatin <- rep(NaN,sum(G))
                         tauhatin[ind0] <- 0
                         fitTauin <- solnp(par = inittauhat, fun=penMSE, b=Btauin[[i]][indnot0], 
                                           A=Ain[[i]][indnot0,indnot0],lam=lambda2,
                                           LB = rep(0,G), eqfun=eqfun, eqB = 1,control=list(trace=0))
                         #tauhatin[indnot0] <- Wminhalf%*%as.vector(fitTauin$pars)
                         tauhatin[indnot0] <- as.vector(fitTauin$pars)
                         
                         return(tauhatin)
                       })
                       tauhatin <- c(list(tauhatin1),tauhatin)
                       
                       ### Compute MSE on left-out part
                       Aouttauin <- lapply(1:nsplits,function(split){Aout[[split]][indnot0,]%*%tauhatin[[split]]})
                       RSStau <- sum(sapply(1:nsplits,function(split){sum((Aouttauin[[split]]-Btauout[[split]][indnot0])^2)/nsplits}))
                       return(RSStau)
                     }
                   }
            )

            #find optimal lambda_2 given muhat
            #first find range for lambda
            minTau <- taus(10^-9) #minimally penalised tau
            maxTau <- taus(10^9) #maximally penalised tau
            lb <- 10^-8
            ub <- 10^8
            diff <- (minTau-maxTau)^2*10^-2 #1 percent relative difference
            while(all(abs(taus(lb)[indnot0]-minTau[indnot0])<diff[indnot0]) & lb<ub){
              lb <- lb*10
            }
            while(all(abs(taus(ub)[indnot0]-maxTau[indnot0])<diff[indnot0]) & ub>lb){
              ub <- ub/10
            }
            rangelambda2 <- c(lb/10,ub*10) #take values just outside the range
            
            #then fit for range of lambda and take minimizer
            if(hypershrinkage=="ridge"){
              lambdas <- 10^seq(log10(rangelambda2[1]),log10(rangelambda2[2]),length.out=100)
            }else{
              lambdas <- 10^seq(log10(rangelambda2[1]),log10(rangelambda2[2]),length.out=30)
            }
            FRSS<-sapply(lambdas,RSSlambdatau)
            minFRSS <- which.min(FRSS)
            if(minFRSS==1) minFRSS <- rev(1:length(lambdas))[which.min(rev(FRSS))] #take least extreme lambda with same RSS
            lambdashat[2] <- lambdas[minFRSS]
            if(profplotRSS){ #profile plot lambda vs RSS
              profPlot <- plot(log10(lambdas),FRSS,xlab="hyperlambda (log10-scale)",ylab="RSS",
                               main=paste("Grouping ",Partitions,", ",hypershrinkage," hypershrinkage",sep=""))
              abline(v=log10(lambdas[minFRSS]),col="red")
              if(!silent) print(paste("Estimated hyperlambda: ",lambdashat[2],sep=""))
            }    
            
            
            # lambda2 <- optimise(RSSlambdatau,rangelambda2) #regular optimiser can get stuck in flat region
            # lambdashat[2] <- lambda2$minimum
            # browser()
          }
        
          #-3.3.3|1.2.4 Compute group variance estimates for optimised hyperpenalty lambda ##############
          if(length(ExtraShrinkage2)==0){
            if(!silent) print(paste("Estimate group weights of grouping ",Partitions,sep=""))
          }
          if(lambdashat[2]==0){
            tauhatold[indnot0] <- solve(A[indnot0,indnot0],Btau[indnot0])
            tauhat <- pmax(0,tauhatold) #set negative tau to 0
          }else{
            tauhatold <- taus(lambdashat[2])
            tauhat <- pmax(0,tauhatold)
            
            if(length(ExtraShrinkage2)>0){
              if(!silent) print(paste("Select groups of grouping ",Partitions,sep=""))
              if(all(tauhatold==0)){ #none selected
                tauhat <- rep(0,G)
                tauhat[indnot0] <- 1
              }else if(sum(tauhatold[indnot0]!=0)==1){ #just one selected
                tauhat <- tauhatold
                Cnorm <- p/sum(c(tauhat)%*%Zt)
                tauhat <- tauhat*Cnorm
              }else{
                #2.
                output <- MoM(Partitions,hypershrinkage=ExtraShrinkage2,
                              pars=list(indnot0=which(tauhatold!=0),ind0=which(tauhatold==0)))
                return(output)
              }
            }
  
          }
          if(normalise){
            Cnorm <- p/sum(c(tauhat)%*%Zt)
            tauhatold <- tauhatold*Cnorm
            weightsTau <- tauhat*Cnorm
            tauhat<-tauhat*Cnorm * tautrgt
          }else{
            weightsTau <- tauhat
            tauhat<- tauhat * tautrgt
            tauhatold <- tauhatold
          }
          
          if(any(is.nan(tauhat))){warning("NaN in group variance");browser()}
        }
      }else{ 
        #-3.3.3|2 Without extra shrinkage --------------------------------------------------------------
        lambdashat <- c(0,0) 
        
        #-3.3.3|2.1 EB estimate group means ============================================================
        muhatp <-as.vector(rep(mu,sum(G))%*%Zt) #px1 vector with estimated prior mean for beta_k, k=1,..,p
        weightsMu <- rep(NaN,sum(G))
        if(!is.nan(mu)){
          muhat<-rep(mu,length(muhat))
        }else{
          if(all(is.nan(betaold))){
            betaold<-rep(1,p) #used as weights
          }else{
            #normalise=F #make sure tau not scaled back to target
          }
          Gamma <- matrix(unlist(
            lapply(Partitions,function(i){ #for each partition
              sapply(1:length(Kg[[i]]),function(j){ #for each group
                #compute row with gamma_{xy}
                x<-groupings[[i]][[j]]
                unlist(sapply(Partitions,function(prt){sapply(groupings[[prt]],function(y){sum(L[x,]%*%t(t(R[,y])/Ik[[prt]][y])%*%diag(betaold[y]))/Kg[[i]][j]})}))
              }, simplify="array")
            })
          ),c(sum(G),sum(G)),byrow=T) #reshape to matrix of size sum(G)xsum(G)
          Bmu <- unlist(
            lapply(Partitions,function(i){ #for each partition
              sapply(1:length(Kg[[i]]),function(j){ #for each group
                x<-groupings[[i]][[j]]
                sum(betasinit[x]-muinitp[x]+L[x,]%*%(R[,pen]%*%muinitp[pen]))/Kg[[i]][j]
              })
            })
          )  
          if(any(is.nan(fixWeightsMu))){ #compute group means for specific partition
            #correct for fixed group means corresponding to groups with variance 0
            if(length(ind0)>0){
              muhat[indnot0] <- solve(Gamma[indnot0,indnot0],Bmu[indnot0]- 
                                        as.matrix(Gamma[indnot0,ind0],c(length(indnot0),length(ind0)))%*%muhat[ind0])
            }else{
              muhat[indnot0] <- solve(Gamma[indnot0,indnot0],Bmu[indnot0])
            }
            #muhat[ind0,Itr+1] <- muhat[ind0,Itr] #means of groups with variance 0 stay the same
            muhatp <-as.vector(c(muhat)%*%Zt)*betaold #px1 vector with estimated prior mean for beta_k, k=1,..,p
            weightsMu <- muhat*p/sum(as.vector(c(muhat)%*%Zt))
          }else{ #compute partition weights/co-data weights
            weightsPart <- sqrt(G[indGrpsGlobal[Partitions]])
            weightMatrixMu <- matrix(rep(0,sum(G)*length(G)),sum(G),length(G))
            for(i in 1:length(G)){
              weightMatrixMu[indGrpsGlobal[[Partitions[i]]],i] <- fixWeightsMu[indGrpsGlobal[[Partitions[i]]]]
            }
            if(!all(round(fixWeightsMu,10)==1)){ #all partitions shrunk to overall mu
              weightsMu <- rep(1/length(Partitions),length(Partitions)) #partition/co-data weights
              muhat<-weightMatrixMu%*%weightsMu #group weights multiplied with partition/co-data weights
            }else{
              Gammatilde <- Gamma%*%weightMatrixMu%*%diag(weightsPart)
              muhat <- solve(t(Gammatilde)%*%Gammatilde,t(Gammatilde)%*%c(Bmu)) / weightsPart
              muhat<- pmax(0,muhat)
              weightsMu <- muhat/sum(muhat) #partition/co-data weights
              muhat<-weightMatrixMu%*%weightsMu #group weights multiplied with partition/co-data weights
            }
          }
          
          
          # if(normalise){ #T by default
          #   C<-mutrgt*p/sum(muhatp)
          #   muhat[,Itr+1]<-muhat[,Itr+1]*C
          #   muhatp <-as.vector(c(muhat[,Itr+1])%*%Zt)*betaold #px1 vector with estimated prior mean for beta_k, k=1,..,p
          # }
          ## Should be same as:
          # Gamma <- ind %*% C %*% diag(betaold) %*% t(ind) /c(Kg)
          # Bmu <- ind %*% (betasinit - Cacc%*%rep(muinit,p)) /c(Kg) #betasinit depend on initial mutrgt
          # muhat <- solve(Gamma,Bmu)
        }
        
        #-3.3.3|2.2 EB estimate group variances ========================================================
        weightsTau <- rep(1,sum(G))
        if(!is.nan(tausq)){
          tauhat <- rep(tautrgt,length(tauhat))
        }else{
          # Btau <- unlist(
          #   lapply(Partitions,function(i){ #for each partition
          #     sapply(1:length(Kg[[i]]),function(j){ #for each group
          #       #compute row with gamma_{xy}
          #       x<-groupings[[i]][[j]]
          #       sum(pmax(0,(betasinit[x]^2-(muinitp[x]+L[x,]%*%(R[,pen]%*%(muhatp[pen]-muinitp[pen])))^2)-V[x]))/Kg[[i]][j]
          #     })
          #   })
          # )
          # A <- matrix(unlist(
          #   lapply(Partitions,function(i){ #for each partition
          #     sapply(1:length(Kg[[i]]),function(j){ #for each group
          #       #compute row with gamma_{xy}
          #       x<-groupings[[i]][[j]]
          #       #compute row with gamma_{xy}
          #       unlist(sapply(Partitions,function(prt){sapply(groupings[[prt]],function(y){
          #         sum(t(L[x,])%*%L[x,]*(R[,y]%*%(t(R[,y])/c(Ik[[prt]][y]))))/Kg[[i]][j]
          #       })}))
          #     }, simplify="array")
          #   })
          # ),c(sum(G),sum(G)),byrow=T) #reshape to matrix of size sum(G)xsum(G)
          # Linear system with division by Variance
          Btau <- unlist(
            lapply(Partitions,function(i){ #for each partition
              sapply(1:length(Kg[[i]]),function(j){ #for each group
                #compute row with gamma_{xy}
                x<-groupings[[i]][[j]]
                sum(pmax(0,(betasinit[x]^2-(muinitp[x]+L[x,]%*%(R[,pen]%*%(muhatp[pen]-muinitp[pen])))^2)/V[x]-1))/Kg[[i]][j]
                #sum((betasinit[x]^2-(muinitp[x]+L[x,]%*%(R[,pen]%*%(muhatp[pen]-muinitp[pen])))^2)/V[x]-1)/Kg[[i]][j]
              })
            })
          )
          A <- matrix(unlist(
            lapply(Partitions,function(i){ #for each partition
              sapply(1:length(Kg[[i]]),function(j){ #for each group
                #compute row with gamma_{xy}
                x<-groupings[[i]][[j]]
                #compute row with gamma_{xy}
                unlist(sapply(Partitions,function(prt){sapply(groupings[[prt]],function(y){
                  sum(t(c(1/V[x])*L[x,])%*%L[x,]*(R[,y]%*%(t(R[,y])/c(Ik[[prt]][y]))))/Kg[[i]][j]
                })}))
              }, simplify="array")
            })
          ),c(sum(G),sum(G)),byrow=T)*tautrgt #reshape to matrix of size sum(G)xsum(G)
          
          if(any(is.nan(fixWeightsTau))){
            if(grepl("positive",hypershrinkage)){
              penMSE <- function(gamma,b,A,lam) return(sum((b-A%*%gamma)^2)) 
              #Aacc <- A%*%Wminhalf
              
              tauhat <- rep(0,G)
              fitTau <- solnp(par = rep(1,length(indnot0)), fun=penMSE, b=Btau[indnot0],
                              A=A[indnot0,indnot0],
                              LB = rep(0,length(indnot0)),control=list(trace=0))
              tauhat[indnot0] <- as.vector(fitTau$pars)

              weightsTau <- tauhat
              tauhat <- tauhat * tautrgt
              tauhatold <- tauhat
            }else{
              tauhat <- rep(0,G)
              tauhatold <- solve(t(A[indnot0,indnot0])%*%A[indnot0,indnot0],
                                 t(A[indnot0,indnot0])%*%Btau[indnot0])
              tauhat <- pmax(0,tauhatold)
              if(normalise){
                #Cnorm <- tautrgt*p/sum(c(tauhat)%*%Zt)
                Cnorm <- p/sum(c(tauhat)%*%Zt)
                weightsTau <- weightsTau*Cnorm
                tauhat<-tauhat*Cnorm * tautrgt
              }else{
                weightsTau <- tauhat
                tauhat<-tauhat * tautrgt
              }
            }
            if(any(is.nan(tauhat))){warning("NaN in group variance")}
          }else{ #compute partition weights/co-data weights
            if(!silent) print("Estimate grouping weights")
            weightsPart <- sqrt(G[Partitions])
            weightMatrixTau <- matrix(rep(0,sum(G)*length(G)),sum(G),length(G))
            for(i in 1:length(G)){
              weightMatrixTau[indGrpsGlobal[[Partitions[i]]],i] <- fixWeightsTau[indGrpsGlobal[[Partitions[i]]]]
            }
            if(all(round(fixWeightsTau,10)==1)){ #all partitions shrunk to overall mu
              weightsTau <- rep(1/length(Partitions),length(Partitions)) #partition/co-data weights
              tauhat<-weightMatrixTau%*%weightsTau*tautrgt #group weights multiplied with partition/co-data weights and overall tau
            }else{
              if(any(partWeightsTau[,Itr]==0)){
                set0 <- unlist(indGrpsGlobal[which(partWeightsTau[,Itr]==0)])
                ind0 <- union(ind0,set0)
                indnot0 <- setdiff(indnot0,set0)
              }
              Atilde <- A[indnot0,indnot0]%*%weightMatrixTau[indnot0,partWeightsTau[,Itr]!=0] 
              
              #Three options to solve for partition weights (use only one):
              #Solve for tau and truncate negative values to 0
              if(1){
                tauhatold<-rep(0,m)
                tauhatold[partWeightsTau[,Itr]!=0] <- solve(t(Atilde)%*%Atilde,t(Atilde)%*%c(Btau[indnot0]))
                weightsTau <- pmax(0,tauhatold)
                #temp<-optim(tauhat,function(x){sum((Atilde%*%x-c(Btau[indnot0]))^2)},lower=rep(0,length(tauhat)),method="L-BFGS-B")
                #tauhat<-temp$par 
                #weightsTau <- tauhat/sum(tauhat) #partition/co-data weights
              }
              
              #solve for w\in[0,1] with convex optimisation package
              #library(CVXR)
              if(0){
                D<-length(G)
                w <- Variable(D)
                objective <- Minimize(sum((Atilde%*%w-c(Btau[indnot0]))^2))
                constraint1 <- diag(rep(1,D))%*%w >= 0
                constraint2 <-  matrix(rep(1,D), nrow = 1)%*%w ==1
                problem <- Problem(objective,constraints = list(constraint1, constraint2))
                result <- solve(problem)
                tauhat <- c(result$getValue(w))
                weightsTau <- pmax(0,tauhat) #correct round-off errors: partition/co-data weights
              }
              
              #Solve for tau>=0 with convex optimisation package, and normalise to 1
              if(0){
                D<-length(G)
                w <- Variable(D)
                objective <- Minimize(sum((Atilde%*%w-c(Btau[indnot0]))^2))
                constraint1 <- diag(rep(1,D))%*%w >= 0
                problem <- Problem(objective,constraints = list(constraint1))
                result <- solve(problem)
                tauhat <- c(result$getValue(w))
                weightsTau <- tauhat#/sum(tauhat) #normalise to 1 to get partition/co-data weights
              }
              
              if(all(tauhat==0)){weightsTau<-rep(0,length(Partitions))}
              if(normalise){
                tauhatold <- tauhatold/sum(weightsTau)
                weightsTau <- weightsTau/sum(weightsTau)
              }
              tauhat<-weightMatrixTau%*%weightsTau*tautrgt #group weights multiplied with partition/co-data weights
            }
          }
        }
      }
      
      return(list(
        lambdashat=lambdashat,
        muhat=muhat,
        tauhat=tauhat,
        tauhatold=tauhatold,
        hypershrinkage=hypershrinkage,
        weightsTau=weightsTau,
        weightsMu=weightsMu
      ))
    }
    
    #For each partition/dataset, use MoM to get group weights
    #NOTE: possible to use different penalty functions
    MoMGroupRes <- lapply(1:m,function(i){
      if(partWeightsTau[i,Itr]!=0){
        MoM(Partitions=i,hypershrinkage=hypershrinkage[i],groupings.grouplvl=groupings.grouplvl[[i]])
      }else{
        return(list(muhat = muhat[indGrpsGlobal[[i]],Itr],
                    tauhatold = tauhatold[indGrpsGlobal[[i]],Itr],
                    tauhat = tauhat[indGrpsGlobal[[i]],Itr],
                    weightsMu = weightsMu[indGrpsGlobal[[i]],Itr],
                    weightsTau = weightsTau[indGrpsGlobal[[i]],Itr],
                    lambdashat = lambdashat[i, Itr,],
                    hypershrinkage=hypershrinkage[i]))
      }
    }
    )

    #global update group parameters
    muhat[,Itr+1]<-unlist(lapply(MoMGroupRes,function(prt){prt$muhat}))
    tauhatold[,Itr+1] <- unlist(lapply(MoMGroupRes,function(prt){prt$tauhatold}))
    tauhat[,Itr+1] <- unlist(lapply(MoMGroupRes,function(prt){prt$tauhat})) #note: scaled version tau
    weightsMu[,Itr+1] <- unlist(lapply(MoMGroupRes,function(prt){prt$weightsMu}))
    weightsTau[,Itr+1] <- unlist(lapply(MoMGroupRes,function(prt){prt$weightsTau}))
    lambdashat[, Itr+1,] <- array(unlist(lapply(MoMGroupRes,function(prt){prt$lambdashat})),c(2,1,m))

    
    #For fixed group weights, use MoM to get partition/co-data weights
    if(m>1){
      if(!is.nan(w)){
        if(is.nan(mu)){
          partWeightsMu[,Itr+1] <- w
          partWeightsMuG[,Itr+1] <- unlist(sapply(1:m,function(x){rep(partWeightsMu[x,Itr+1],G[x])})) #total number of groups x 1 vector with partition weights
        }
        if(is.nan(tausq)){
          partWeightsTau[,Itr+1] <- w
          partWeightsTauG[,Itr+1] <- unlist(sapply(1:m,function(x){rep(partWeightsTau[x,Itr+1],G[x])})) #total number of groups x 1 vector with partition weights
        }
      }else{
        if(!all(round(weightsTau[,Itr+1],10)==1)){
          #if(!any(partWeightsTau[,Itr]==0)){
            MoMPartRes <- MoM(Partitions=1:m,hypershrinkage="none",fixWeightsMu=weightsMu[,Itr+1],fixWeightsTau=weightsTau[,Itr+1])
          # }else{
          #   partNot0 <- which(partWeightsTau[,Itr]!=0)
          #   MoMPartRes <- MoM(Partitions=partNot0,hypershrinkage="none",
          #                     fixWeightsMu=weightsMu[unlist(indGrpsGlobal[partNot0]),Itr+1],
          #                     fixWeightsTau=weightsTau[unlist(indGrpsGlobal[partNot0]),Itr+1])
          # }
          
        }
        if(is.nan(mu)){
          partWeightsMu[,Itr+1] <- MoMPartRes$weightsMu
          partWeightsMuG[,Itr+1] <- unlist(sapply(1:m,function(x){rep(partWeightsMu[x,Itr+1],G[x])})) #total number of groups x 1 vector with partition weights
        }
        if(is.nan(tausq)){
          if(all(round(weightsTau[,Itr+1],10)==1)){
            partWeightsTau[,Itr+1] <- rep(1/m,m)
            partWeightsTauG[,Itr+1] <- unlist(sapply(1:m,function(x){rep(partWeightsTau[x,Itr+1],G[x])})) #total number of groups x 1 vector with partition weights
          }else{
            partWeightsTau[,Itr+1] <- pmax(MoMPartRes$weightsTau,0)
            partWeightsTauG[,Itr+1] <- unlist(sapply(1:m,function(x){rep(partWeightsTau[x,Itr+1],G[x])})) #total number of groups x 1 vector with partition weights
          }
        }
      }
    }

    if(all(is.nan(betaold))){
      betaold <-rep(1,p) #px1 vector with estimated prior mean for beta_k, k=1,..,p
    }
    if(is.nan(mu)){
      muhatp <- as.vector(c(partWeightsMuG[,Itr+1]*muhat[,Itr+1])%*%Zt)*betaold[pen]
    }else{
      muhatp<-rep(mu,length(pen))
    }
    
    #-3.3.4 Update group-specific penalties ###################################################################
    if(is.nan(tausq)){
      if(all(partWeightsTauG==0)){#set all partition/group weights to 1 (i.e. no difference in partitions/groups)
        lambdap<-sigmahat/as.vector(c(tauhat[,1])%*%Zt) #target tau/overall
        }else{
          lambdap<-sigmahat/as.vector(c(partWeightsTauG[,Itr+1]*tauhat[,Itr+1])%*%Zt) #specific penalty for beta_k
          lambdap[lambdap<0]<-Inf 
      } 
    }else{
      lambdap<-rep(sigmahat/tausq,length(pen))
    }
    #should be the same as
    #lambdap2<-sigmahat/as.vector(c(partWeightsTauG*weightsTau*tautrgt)%*%Zt)
  
    #-3.3.5 Update beta using glmnet #######################################################################
    if(!silent) print("Estimate regression coefficients")
    if(all(tauhat[,Itr+1]==0)){
      beta <- muhatp
      if(intrcpt){
        if(model=="linear"){
          glmGR <- list(a0=sum(Y-X%*%beta)/n)
        }else if(model=='logistic'){
          glmGR <- list(a0=sum(Y-exp(X%*%beta)/(1+exp(X%*%beta)))/n) 
        }
      }else{
        glmGR <- list(a0=0)
      }
      warning("All tau (set to) 0")
    }else{
      if(model=="linear"){
        sd_y2 <- sqrt(var(Y-X %*% muhatp)*(n-1)/n)[1]
      }else if(model%in%c('logistic','cox')){
        sd_y2 <- 1 #do not standardize Y-offset for logistic/cox model
      }
      if(any(is.nan(sqrt(lambdap)))){browser()}
      lambdaoverall <- sigmahat/tautrgt
      Xacc <- X
      Xacc[,pen] <- as.matrix(X[,pen] %*% sparseMatrix(i=1:length(pen),j=1:length(pen),
                                                       x=c(1/sqrt(lambdap/lambdaoverall))))
      if(model=="cox"){
        glmGR <- glmnet(Xacc,Y,alpha=0,
                        lambda = lambdaoverall/n*sd_y2,family=fml,
                        offset = X[,!((1:p)%in%unpen)] %*% muhatp, standardize = F,
                        penalty.factor=penfctr, thresh=10^-10)
        # glmGR <- glmnet(X,Y,alpha=0,
        #                 lambda = lambdaoverall/n*sd_y2,family=fml,
        #                 offset = X[,!((1:p)%in%unpen)] %*% muhatp, standardize = F,
        #                 penalty.factor=penfctr*lambdap/lambdaoverall, thresh=10^-10)
      }else{
        glmGR <- glmnet(Xacc,Y,alpha=0,
                        lambda = lambdaoverall/n*sd_y2,family=fml,
                        offset = X[,!((1:p)%in%unpen)] %*% muhatp, intercept = intrcpt, standardize = F,
                        penalty.factor=penfctr, thresh=10^-10)
        # glmGR <- glmnet(X,Y,alpha=0,
        #                 lambda = lambdaoverall/n*sd_y2,family=fml,
        #                 offset = X[,!((1:p)%in%unpen)] %*% muhatp, intercept = intrcpt, standardize = F,
        #                 penalty.factor=penfctr*lambdap/lambdaoverall, thresh=10^-10)
      }
      beta <- as.vector(glmGR$beta) 
      beta[pen] <- c(1/sqrt(lambdap/lambdaoverall)) * beta[pen] + muhatp
      # beta[pen] <- as.vector(glmGR$beta)[pen] + muhatp
    }
    
    #-3.3.6 Update predictions on independent data (if given) ################################################
    if(!missing(X2)){
      #Ypredridge <- predict(glmGR,newx=X2)
      if(model=="linear"){
        X2c <- cbind(X2,rep(1,n2))
        YpredGR[,Itr+1] <- X2c %*% c(beta,glmGR$a0)
        MSEecpc[Itr+1]<- sum((YpredGR[,Itr+1]-Y2)^2)/n2
      } 
      if(model=='logistic'){
        X2c <- cbind(X2,rep(1,n2))
        YpredGR[,Itr+1] <- 1/(1+exp(-X2c %*% c(beta,glmGR$a0)))
        MSEecpc[Itr+1]<- sum((YpredGR[,Itr+1]-Y2)^2)/n2
        if(any(is.nan(YpredGR[,Itr+1]))){browser()}
      }else if(model=='cox'){ 
        expXb<-exp(X %*% c(beta))
        h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXb[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
        H0 <- sapply(Y2[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
        YpredGR[,Itr+1] <- H0*exp(X2 %*% beta)
        MSEecpc[Itr+1]<- sum((YpredGR[,Itr+1]-Y2[,2])^2)/n2
      }
    }
    
    #-3.3.7 Set current estimates as initial estimates for next iteration #####################################
    betasinit <- beta
    muinitp <- muhatp
    intrcptinit <- glmGR$a0
    ind0 <- which(tauhat[,Itr+1]==0) #update index of groups with zero variance
    indnot0 <- which(tauhat[,Itr+1]>0)
    Itr <- Itr+1

    if(length(ind0)==sum(G)){
      if(!silent) print(paste("All prior group variances estimated to be 0, iterating stopped after",Itr,"iteration(s)"))
      break;
    }
  }
  
  #-4. (optional) Posterior variable selection-------------------------------------------------------
  #postselection: indicates which of the two methods possible is used to do post-hoc variable selection
  #maxsel: maximum number of variables selected
  if(postselection!=F){
    if(!silent) print("Sparsify model with posterior selection")
    postSel <- postSelect(X=X,Y=Y,beta=beta,intrcpt=glmGR$a0,penfctr=penfctr, 
                             postselection=postselection,maxsel=maxsel, 
                             penalties=lambdap,model=model,tauglobal=tautrgt,
                             sigmahat=sigmahat,muhatp=muhatp, 
                             X2=X2,Y2=Y2)
  }
  
  #-5. (optional) Compare with glmnet --------------------------------------------------------------
  betaridge<-NaN; Ypredridge<-NaN; 
  if(!is.nan(compare) & compare){
    if(model=="cox"){
      glmR <- glmnet(X,Y,family=fml,alpha=0,
                     lambda=lambdaridge*sd_y/n,standardize = F,
                     penalty.factor=penfctr)
    }else{
      glmR <- glmnet(X,Y,family=fml,alpha=0,
                     lambda=lambdaridge*sd_y/n,standardize = F,intercept=intrcptGLM,
                     penalty.factor=penfctr)
    }
    betaridge <- as.vector(glmR$beta)
    
    if(!missing(X2)){
      #Ypredridge <- predict(glmR,newx=X2)
      #browser()
      
      if(model=="linear"){
        X2c <- cbind(X2,rep(1,n2))
        Ypredridge <- X2c %*% c(betaridge,glmR$a0)
        MSEridge <- sum((Ypredridge-Y2)^2)/n2
      } 
      if(model=='logistic'){
        X2c <- cbind(X2,rep(1,n2))
        Ypredridge <- 1/(1+exp(-X2c %*% c(betaridge,glmR$a0)))
        MSEridge <- sum((Ypredridge-Y2)^2)/n2
      }else if(model=="cox"){
        expXb<-exp(X %*% c(betaridge))
        h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXb[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
        H0 <- sapply(Y2[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
        Ypredridge <- H0*exp(X2 %*% betaridge)
        MSEridge<- sum((Ypredridge-Y2[,2])^2)/n2
      }
    }
  }

  #-6. Output -------------------------------------------------------------------------------------
  output <- list(
    beta=beta, #beta from ecpc (with Group Ridge penalties)
    intercept=intrcptinit, #unpenalised intercept covariate
    tauglobal=tautrgt, #overall tautrgt
    gammatilde = tauhatold[,nIt+1], #EB estimated prior group variance before truncating
    gamma=weightsTau[,nIt+1], #group weights variance
    w = partWeightsTau[,nIt+1], #grouping weights in local variances
    penalties = lambdap, #penalty parameter on all p covariates
    hyperlambdas = lambdashat[2,nIt+1,], #hyperpenalties for all groupings
    #weights = weights, #weights used in ridge hypershrinkage
    #levelsY = levelsY, #in case of logistic
    sigmahat=sigmahat #estimated sigma^2 (linear model)
  )
  if(nIt>1){
    output$gamma <- weightsTau[,-1]
    output$gammatilde <- tauhatold[,-1]
    output$w <- partWeightsTau[,-1]
    output$hyperlambdas <- lambdashat[2,-1,]
  }
  if(!missing(X2)){
    output$Ypred<-YpredGR[,-1] #predictions for test set
    output$MSEecpc <- MSEecpc[nIt+1] #MSE on test set
    if(nIt>1){
      output$Ypred<-YpredGR[,-1] #predictions for test set
      output$MSEecpc <- MSEecpc[-1] #MSE on test set
    }
  }
  
  if(mutrgt!=0){ #prior group means are estimated as well
    output$muhat <- muhat[,nIt+1] #EB estimated prior group means: mu_global*muweight
    output$muglobal <- mutrgt #overall mutrgt
    output$gamma.mu <- weightsMu[,nIt+1] #group weights mean
    output$w.mu <- partWeightsMu[,nIt+1] #grouping weights mean
    output$hyperlambdas.mu <- lambdashat[1,nIt+1,] #hyperpenalties for all groupings for the group means
    output$offset.lp <- X[,!((1:p)%in%unpen)] %*% muhatp #offset used in computing final beta (default 0),
  }

  if(!is.nan(compare) & compare){ #comparison with ordinary ridge obtained with glmnet
    output$betaridge <- betaridge #ordinary ridge beta
    output$interceptridge <- glmR$a0
    output$lambdaridge <- lambdaridge #ordinary ridge lambda
    if(!all(is.nan(X2))){
      output$Ypredridge <- Ypredridge
      output$MSEridge <- MSEridge
    }
  }
  if(postselection!=F){ #posterior selection is performed
    output$betaPost <- postSel$betaPost
    output$interceptPost <- postSel$a0
    if(!missing(X2)){
      output$MSEPost <- postSel$MSEPost #MSE on independent data set (if given)
      output$YpredPost <- postSel$YpredPost #predictions for independent data set (if given)
    }
  }
  return(output)
}

### Other functions 
#Select covariates a posteriori----
postSelect <- function(X,Y,beta,intrcpt=0,penfctr, #input data
                          postselection="elnet+dense",maxsel=30, #selection options
                          penalties,model,tauglobal,sigmahat=1,muhatp=0, #needed for method "elnet"
                          X2=NaN,Y2=NaN){
  #Description:
  #Post-hoc variable selection: select maximum maxsel covariates of all penalised covariates
  #Unpenalised covariates (e.g. intercept) are always selected, on top of the maxsel number of selected penalised covariates
  #
  #Input data:
  #-X: (nxp) data matrix: data of p penalised and unpenalised covariates on n samples
  #-Y: (nx1) vector: response
  #-beta: previous estimate in dense model
  #-intrcpt: value of intercept in dense model
  #-penfctr: as in glmnet penalty.factor. Default: 0 if covariate is not penalised, 1 if covariate is penalised
  #
  #Input options:
  #-postselection: "elnet,dense" (default), "elnet,sparse", "BRmarginal,dense", "BRmarginal,sparse" or "DSS" 
  #               for corresponding post-selection method used and dense/sparse indicating if the global tau^2 is recalibrated 
  #              (for sparse settings) or not (for dense settings)
  #-maxsel: maximum number of penalised covariates to be selected (additional to all unpenalised covariates)
  #
  #Input additional parameters method "elnet":
  #-penalties: (length(pen)x1) vector: ridge penalties for all penalised covariates found in dense model
  #-model: "linear", "logistic" or "cox" for type of regression model used
  #-tauglobal, sigmahat, muhatp: parameter estimates of dense model fit
  #
  #Input optional:
  #X2,Y2 (optional): independent data and response on which predictions and MSE is computed
  tautrgt <- tauglobal 
  n<-dim(X)[1] #number of samples
  p<-dim(X)[2] #number of covariates (penalised and unpenalised)
  if(missing(penfctr)) penfctr <- rep(1,p) #all covariates penalised the same

  maxsel2 <- pmin(maxsel, p)
  if(any(maxsel2<2)){
    warning("Number of variables to be selected should be at least 2 (out of convenience)")
    if(!silent) print("Maxsel values smaller than 2 are set to 2")
    maxsel2[maxsel2<2] <- 2
  }
  nonzeros <- beta!=0 #Fix beta that are already 0 
  lambdap<-penalties #ridge penalties for the penalised and non-zero covariates
  pen<-which(penfctr!=0) #index of penalised covariates

  switch(model,
         'linear'={
           fml <- 'gaussian'
           sd_y <- sqrt(var(Y)*(n-1)/n)[1]
         },
         'logistic'={
           fml <- 'binomial'
           sd_y <- 1 #do not standardise y in logistic setting
         },
         'cox'={
           fml <- 'cox'
           sd_y <- 1 #do not standardise y in cox regression setting
         }
  )
  
  if(grepl("elnet",postselection)){
    #multiple numbers possible of maximum number of covariates to be selected
    output<-lapply(maxsel2,function(x){
      
      #if already less than maxsel covariates selected
      if(sum(beta[pen]!=0)<=x){
        betaPost <- beta
        
        output<-list()
        output$betaPost <- betaPost #post-selection beta using Elastic Net
        output$whichPost <- which(nonzeros & (1:p)%in%pen) #index of selected penalised covariates
        output$a0 <- intrcpt
        
        if(!all(is.nan(X2))){
          if(model=="linear"){
            YpredPost <- X2 %*% c(betaPost) + intrcpt
            MSEPost <- sum((YpredPost-Y2)^2)/length(Y2)
          } 
          if(model=='logistic'){
            YpredPost <- 1/(1+exp(-(X2 %*% c(betaPost) + intrcpt)))
            MSEPost <- sum((YpredPost-Y2)^2)/length(Y2)
            if(any(is.nan(YpredPost))){browser()}
          }else if(model=='cox'){ 
            expXbPost<-exp(X %*% c(betaPost))
            h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXbPost[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
            H0 <- sapply(Y2[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
            YpredPost <- H0*exp(X2 %*% c(betaPost))
            MSEPost<- sum((YpredPost-Y2[,2])^2)/length(Y2[,2])
          }
          
          output$MSEPost <- MSEPost #MSE on independent data set (if given)
          output$YpredPost <- YpredPost #predictions for independent data set (if given)
        }
        
        return(output)
      }else{
        if(length(intrcpt)==0||intrcpt==0){intrcpt <- F}else{intrcpt <- T}
        if(length(muhatp)==1) muhatp <- rep(muhatp,length(pen))
        if(all(muhatp==0)) offset <- rep(0,n)
        else offset <- X[,pen] %*% muhatp #in case prior mean of penalised covariates is not equal to 0
        
        lambdaoverall <- sigmahat/tautrgt
        lam2 <- sigmahat/tautrgt/n*sd_y
        Xacc <- X
        Xacc[,pen] <- as.matrix(X[,pen] %*% sparseMatrix(i=1:length(lambdap),j=1:length(lambdap),
                                                   x=c(1/sqrt(lambdap/lambdaoverall))) )

        #define function with output number of selected variables minus maximum possible
        #find root of function such that we have at most maxsel variables
        fsel <- function(alpha, maxselec = x) {
          if(alpha == 0) 
            return(p - maxselec)
          else {
            if(model=="cox"){
              glmPost <- glmnet(Xacc[,nonzeros],Y,alpha=alpha,
                              lambda = lam2/(1-alpha),family=fml,
                              offset = offset, standardize = F,
                              penalty.factor=penfctr[nonzeros], thresh=10^-10)
            }else if(model %in% c("logistic","linear")){
              glmPost <- glmnet(Xacc[,nonzeros],Y,alpha=alpha,
                              lambda = lam2/(1-alpha),family=fml,
                              offset = offset, standardize = F,
                              intercept= intrcpt,
                              penalty.factor=penfctr[nonzeros], thresh=10^-10)
            }
            betaPost <- rep(0,p)
            betaPost[nonzeros] <- as.vector(glmPost$beta)
            betaPost[pen] <- c(1/sqrt(lambdap/lambdaoverall)) * betaPost[pen] + muhatp 
            return(sum(betaPost[pen] != 0) - maxselec ) #number of non-zero penalised covariates
          }
        }
        #Elastic net uses both lasso penalty lam1 and ridge penalty lam2
        #Keep lam2 fixed, use alpha and lambda arguments of glmnet to search in range lam1\in[0,10*lam2]
        #Equivalent to searching in alpha\in[0,10/11] and corresponding lambda
        rangeAlpha <- c(0,10/11)
        ItrAlp <- 1
        #if (sign(fsel(0))==sign(fsel(10/11))) browser()
        while(sign(fsel(0))==sign(fsel(rangeAlpha[2])) & ItrAlp <=50){
          rangeAlpha[2] <- rangeAlpha[2] + 0.5*(1-rangeAlpha[2])
          ItrAlp <- ItrAlp +1 
        }
        alpha <- uniroot(fsel, interval = rangeAlpha, maxiter = 200,tol=10^(-10))$root
        
        #for found alpha, refit model to see which beta are selected
        if(model=="cox"){
          glmPost0 <- glmnet(Xacc[,nonzeros],Y,alpha=alpha,
                           lambda = lam2/(1-alpha),family=fml,
                           offset = offset, standardize = F,
                           penalty.factor=penfctr[nonzeros], thresh=10^-10)
        }else{
          glmPost0 <- glmnet(Xacc[,nonzeros],Y,alpha=alpha,
                           lambda = lam2/(1-alpha) ,family=fml,
                           offset = offset, intercept = intrcpt, standardize = F,
                           penalty.factor=penfctr[nonzeros], thresh=10^-10)
        }
        betaPost0 <- rep(0,p)
        betaPost0[nonzeros] <- as.vector(glmPost0$beta)
        betaPost0[pen] <- c(1/sqrt(lambdap/sigmahat*tautrgt)) * betaPost0[pen] + muhatp 
        whichPostboth <- betaPost0 != 0 #both unpenalised covariates and selected penalised

        if(sum(whichPostboth)<=1){
          warning("At least two variables should be selected for glmnet")
          return(list(betaPost=NULL,whichPost=NULL,a0=NULL))
        }
        if(grepl("dense",postselection)){ #use weighted penalty
          if(!all(muhatp==0)) offset <- X[,whichPostboth & (1:p)%in%pen, drop=F] %*% muhatp[whichPostboth[pen], drop=F]
          
          #recalibrate overall lambda using cross-validation on selected variables only
          if(grepl("dense2",postselection)){ 
            if(model=="cox"){
              lambdaGLM<-cv.glmnet(Xacc[,whichPostboth, drop=F],Y,alpha=0,
                                   family=fml,offset = offset ,standardize = F,
                                   penalty.factor=penfctr[whichPostboth], thresh=10^-10) #alpha=0 for ridge
              lam2<-lambdaGLM$lambda.min
            }else{
              lambdaGLM<-cv.glmnet(Xacc[,whichPostboth, drop=F],Y,alpha=0,
                                   family=fml,offset = offset , intercept = intrcpt, standardize = F,
                                   penalty.factor=penfctr[whichPostboth], thresh=10^-10) #alpha=0 for ridge
              lam2<-lambdaGLM$lambda.min
            }
          }
          
          #Recompute beta using only selected beta and group ridge penalty (without lasso penalty)
          if(model=="cox"){
            glmPost <- glmnet(Xacc[,whichPostboth, drop=F],Y,alpha=0,
                            lambda = lam2,family=fml,
                            offset = offset ,standardize = F,
                            penalty.factor=penfctr[whichPostboth], thresh=10^-10)
          }else{
            glmPost <- glmnet(Xacc[,whichPostboth, drop=F],Y,alpha=0,
                            lambda = lam2,family=fml,
                            offset = offset , intercept = intrcpt, standardize = F,
                            penalty.factor=penfctr[whichPostboth], thresh=10^-10)
          }
          betaPost <- rep(0,p)
          betaPost[whichPostboth] <- as.vector(glmPost$beta)
          indPostpen <- whichPostboth[pen] #penalised covariates in range 1:length(pen) that are selected and non-zero
          indPostp <- whichPostboth&((1:p)%in%pen) #penalised covariates in range 1:p that are selected and non-zero
          betaPost[indPostp] <- c(1/sqrt(lambdap[indPostpen]/sigmahat*tautrgt)) * betaPost[indPostp] + muhatp[indPostpen]
          
          whichPost <- which(indPostp) #index of selected penalised covariates 
          output<-list()
          output$betaPost <- betaPost #post-selection beta using Elastic Net
          output$whichPost <- whichPost #index of selected covariates
          output$a0 <- glmPost$a0
          #output$offsetPost <- offset #offset used in Post
        }else{# if(grepl("sparse",postselection)){ #refit standard ridge with newly cross-validated lambda
          if(model=="cox"){
            lambdaGLM<-cv.glmnet(X[,whichPostboth, drop=F],Y,alpha=0,family=fml,
                                 standardize = F,penalty.factor=penfctr[whichPostboth]) #alpha=0 for ridge
            lam2<-lambdaGLM$lambda.min
          }else{
            lambdaGLM<-cv.glmnet(X[,whichPostboth, drop=F],Y,alpha=0,family=fml,
                                 standardize = F,penalty.factor=penfctr[whichPostboth]) #alpha=0 for ridge
            lam2<-lambdaGLM$lambda.min
          }
          if(model=="cox"){
            glmPost <- glmnet(X[,whichPostboth, drop=F],Y,alpha=0,
                            lambda = lam2,family=fml,
                            offset = offset ,standardize = F,
                            penalty.factor=penfctr[whichPostboth], thresh=10^-10)
          }else{
            glmPost <- glmnet(X[,whichPostboth, drop=F],Y,alpha=0,
                            lambda = lam2,family=fml,
                            intercept = intrcpt, standardize = F,
                            penalty.factor=penfctr[whichPostboth], thresh=10^-10)
          }
          betaPost <- rep(0,p)
          betaPost[whichPostboth] <- as.vector(glmPost$beta)
          indPostpen <- whichPostboth[pen] #penalised covariates in range 1:length(pen) that are selected and non-zero
          indPostp <- whichPostboth&((1:p)%in%pen) #penalised covariates in range 1:p that are selected and non-zero
          
          whichPost <- which(indPostp) #index of selected penalised covariates 
          output<-list()
          output$betaPost <- betaPost #post-selection beta using Elastic Net
          output$whichPost <- whichPost #index of selected covariates
          output$a0 <- glmPost$a0
          #output$offsetPost <- offset #offset used in Post
        }
        
        if(!all(is.nan(X2))){
          if(model=="linear"){
            YpredPost <- X2 %*% c(betaPost) + output$a0
            MSEPost <- sum((YpredPost-Y2)^2)/length(Y2)
          } 
          if(model=='logistic'){
            YpredPost <- 1/(1+exp(-(X2 %*% c(betaPost) + output$a0)))
            MSEPost <- sum((YpredPost-Y2)^2)/length(Y2)
            if(any(is.nan(YpredPost))){browser()}
          }else if(model=='cox'){ 
            expXbPost<-exp(X %*% c(betaPost))
            h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXbPost[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
            H0 <- sapply(Y2[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
            YpredPost <- H0*exp(X2 %*% c(betaPost))
            MSEPost <- sum((YpredPost-Y2[,2])^2)/length(Y2[,2])
          }
          
          output$MSEPost <- MSEPost #MSE on independent data set (if given)
          output$YpredPost <- YpredPost #predictions for independent data set (if given)
        }
        
        return(output)
      }
    })
  }else if(grepl("DSS",postselection)){
      #Hahn and Carvalho 2014
      #Use adaptive lasso penalty on linear predictiors as in Equation (22) of their paper
      Xacc <- X
      Xacc[,pen] <- as.matrix(Xacc[,pen]%*% sparseMatrix(i=1:length(beta[pen]),j=1:length(pen),
                                   x=c(sqrt(abs(beta[pen])))) )
      Ygamma <- Xacc%*%beta
      
      
      if(grepl("fast",postselection)){
        #use glmnet to compute beta for a whole range of lambda simultaneously
        #faster, but might result in fewer selected variables than asked for
        glmPost <- glmnet(x=Xacc[,nonzeros],y=Ygamma,alpha=1,nlambda=100,
                          family="gaussian",
                          intercept = intrcpt, standardize = F,
                          penalty.factor=penfctr[nonzeros], thresh=10^-10)
      }else{ 
        #retrieve exactly maxsel non-zero covariates
        if(length(intrcpt==0)||intrcpt==0){intrcpt <- F}else{intrcpt <- T}
        fsel <- function(lam1,maxselec=maxsel2){
          glmPost <- glmnet(x=Xacc[,nonzeros],y=Ygamma,alpha=1,
                            lambda = lam1,family="gaussian",
                            intercept = intrcpt, standardize = F,
                            penalty.factor=penfctr[nonzeros], thresh=10^-10)
          return(sum(as.vector(glmPost$beta) !=0) - maxselec)
        }
      }
      
      output <- lapply(maxsel2,function(x){
        #if already less than maxsel covariates selected
        if(length(beta[pen]!=0)<=x){
          betaPost <- beta
          
          output<-list()
          output$betaPost <- betaPost #post-selection beta using Elastic Net
          output$whichPost <- which(betaPost!=0) #index of selected covariates
          output$a0 <- intrcpt
          
          if(!all(is.nan(X2))){
            if(model=="linear"){
              YpredPost <- X2 %*% c(betaPost) + intrcpt
              MSEPost <- sum((YpredPost-Y2)^2)/length(Y2)
            } 
            if(model=='logistic'){
              YpredPost <- 1/(1+exp(-(X2 %*% c(betaPost) + intrcpt)))
              MSEPost <- sum((YpredPost-Y2)^2)/length(Y2)
              if(any(is.nan(YpredPost))){browser()}
            }else if(model=='cox'){ 
              expXbPost<-exp(X %*% c(betaPost))
              h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXbPost[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
              H0 <- sapply(Y2[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
              YpredPost <- H0*exp(X2 %*% c(betaPost))
              MSEPost <- sum((YpredPost-Y2[,2])^2)/length(Y2[,2])
            }
            
            output$MSEPost <- MSEPost #MSE on independent data set (if given)
            output$YpredPost <- YpredPost #predictions for independent data set (if given)
          }
          
          return(output)
        }else{
          if(grepl("fast",postselection)){
            whlam <- which(glmPost$df<=x)
            takelam <- rev(whlam)[1]
            lam <- glmPost$lambda[takelam]
            
            betaPost <- rep(0,p)
            betaPost[nonzeros]<- as.vector(glmPost$beta[,takelam])
            betaPost[pen] <- c(sqrt(abs(beta[pen]))) * betaPost[pen]
            whichPost <- which(betaPost!=0 & ((1:p)%in%pen))
            a0 <- glmPost$a0[takelam]
          }else{
            lam1 <- uniroot(fsel,maxselec=x, interval = c(0, max(lambdap[lambdap<Inf])*sigmahat/n), 
                            maxiter = 200,tol=10^(-10))$root
            glmPost <- glmnet(x=Xacc[,nonzeros],y=Ygamma,alpha=1,
                              lambda = lam1,family="gaussian",
                              intercept = intrcpt, standardize = F,
                              penalty.factor=penfctr[nonzeros], thresh=10^-10)
            betaPost <- rep(0,p)
            betaPost[nonzeros]<- as.vector(glmPost$beta)
            betaPost[pen] <- c(sqrt(abs(beta[pen]))) * betaPost[pen]
            whichPost <- which(betaPost!=0 & ((1:p)%in%pen))
            a0<-glmPost$a0
          }
          
          output<-list()
          output$betaPost <- betaPost #post-selection beta using Elastic Net
          output$whichPost <- whichPost #index of selected covariates
          output$a0 <- a0
        
          if(!all(is.nan(X2))){
            if(model=="linear"){
              YpredPost <- X2 %*% c(betaPost) + output$a0
              MSEPost<- sum((YpredPost-Y2)^2)/length(Y2)
            } 
            if(model=='logistic'){
              YpredPost <- 1/(1+exp(-(X2 %*% c(betaPost) + output$a0)))
              MSEPost<- sum((YpredPost-Y2)^2)/length(Y2)
              if(any(is.nan(YpredPost))){browser()}
            }else if(model=='cox'){ 
              expXbPost<-exp(X %*% c(betaPost))
              h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXbPost[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
              H0 <- sapply(Y2[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
              YpredPost <- H0*exp(X2 %*% c(betaPost))
              MSEPost<- sum((YpredPost-Y2[,2])^2)/length(Y2[,2])
            }
            
            output$MSEPost <- MSEPost #MSE on independent data set (if given)
            output$YpredPost <- YpredPost #predictions for independent data set (if given)
          }
          return(output)
        }
    })
  }else if(grepl("BR",postselection)){
    #Bondell & Reich propose two options:
    #1. (joint): joint credible set approach using posterior covariance matrix
    #2. (marginal): marginal credible set approach using posterior standard deviations
    #After selecting, predict with following options:
    #1. (same): Keep output selection procedure
    #2. (dense): Refit with previously estimated weighted ridge
    #3. (sparse): Refit using ordinary ridge with newly cross-validated lambda 
    if(grepl("joint",postselection)){
      #joint credible set approach
      if(!silent) print("Bondell&Reich joint credible set used to select covariates")
      
      ind <- which(penfctr!=0 & penalties!=Inf) #index of beta that are penalised and not already set to 0
      #compute prediction weight matrix W
      if(model=="linear"){
        W<-diag(rep(1,n))
      }else if(model=="logistic"){
        expminXb<-exp(-X%*%c(beta) - intrcpt)
        Pi<-1/(1+expminXb)
        W<-diag(c(sqrt(Pi*(1-Pi))))
      }else if(model=="cox"){
        expXb<-exp(X%*%c(beta))
        h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXb[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
        H0 <- sapply(Y[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
        W <- diag(c(sqrt(H0*expXb)))
      }
      
      Sigminhalf<-sqrtm(t(X[,ind])%*%W%*%X[,ind]+diag(lambdap[ind])) / sqrt(sigmahat) #posterior covariance matrix Sigma, Sigma^{-0.5}
      #D<- diag(beta[ind]^2) #diagonal matrix with beta_k^2 on the diagonal
      
      Xstar <- t(t(Sigminhalf) * beta[ind]^2) #Sigma^{-0.5} D
      Ystar <- Sigminhalf %*% beta[ind] 
      
      #first fit for range of lambda 
      glmPost <- glmnet(Xstar,Ystar,alpha=1,family="gaussian",
                      standardize = F,intercept= F, thresh=10^-10)
     
      output<-lapply(maxsel2,function(x){
        fsel <- function(lambda,maxselec = x){
          fit<-coef.glmnet(glmPost,s=lambda,exact=T,x=Xstar,y=Ystar)
          return(sum(fit!=0)- maxselec) #number of non-zero penalised covariates
        }
        lambda <- uniroot(fsel, interval = c(0, max(glmPost$lambda)), maxiter = 200,tol=10^(-10))$root
        fit<-coef.glmnet(glmPost,s=lambda,exact=T,x=Xstar,y=Ystar)
      
        if(grepl("same",postselection)){
          if(x==maxsel2[1]) if(!silent) print("Selected covariates are not refit, unpenalised covariates are kept the same")
          betaPost <- beta
          betaPost[ind] <- fit*beta[ind]^2
          whichPost <- which(betaPost!=0 & penfctr!=0) #index of selected penalised covariates 
          output<-list()
          output$betaPost <- betaPost #post-selection beta using Elastic Net
          output$whichPost <- whichPost #index of selected covariates
          output$a0 <- intrcpt
          
          if(!all(is.nan(X2))){
            if(model=="linear"){
              YpredPost <- X2 %*% c(betaPost) + intrcpt
              MSEPost <- sum((YpredPost-Y2)^2)/length(Y2)
            } 
            if(model=='logistic'){
              YpredPost <- 1/(1+exp(-(X2 %*% c(betaPost) + intrcpt)))
              MSEPost <- sum((YpredPost-Y2)^2)/length(Y2)
              if(any(is.nan(YpredPost))){browser()}
            }else if(model=='cox'){ 
              expXbPost<-exp(X %*% c(betaPost))
              h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXbPost[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
              H0 <- sapply(Y2[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
              YpredPost <- H0*exp(X2 %*% c(betaPost))
              MSEPost<- sum((YpredPost-Y2[,2])^2)/length(Y2[,2])
            }
          
            output$MSEPost <- MSEPost #MSE on independent data set (if given)
            output$YpredPost <- YpredPost #predictions for independent data set (if given)
          }
          return(output)
          }else if(grepl("sparse",postselection)){
            if(!silent) print("Selected covariates are refit with previously estimated weighted ridge prior")
            
          }else{# if(grepl("dense",postselection)){
            if(!silent) print("Selected covariates are refit with an ordinary ridge prior using newly cross-validated penalty")
            
          }
        })
    }else{ #if(grepl("marginal")){
      #marginal credible set approach
      if(!silent) print("Bondell&Reich marginal credible set used to select covariates")
      
      ind <- which(penfctr!=0 & beta!=0) #index of betas that are penalised and not already set to 0
      indpen <- which((beta[penfctr!=0])!=0) #index in 1:length(pen) of betas that are penalise and not yet set to 0
      #compute prediction weight matrix W
      if(model=="linear"){
        diagWhalf<-rep(1,n)
      }else if(model=="logistic"){
        expminXb<-exp(-X%*%c(beta) - intrcpt)
        Pi<-1/(1+expminXb)
        diagWhalf<-c(Pi*(1-Pi)) #W^0.5
      }else if(model=="cox"){
        expXb<-exp(X%*%c(beta))
        h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXb[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
        H0 <- sapply(Y[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
        diagWhalf <- c(H0*expXb)
      }
      Xtilde <- t(t(diagWhalf * X[,ind]) / sqrt(lambdap[indpen]))
      svdXtilde <- svd(Xtilde)
      sdBetas <- 1/sqrt(lambdap[indpen]) * sqrt(1-diag(svdXtilde$v %*% (1/(1+1/svdXtilde$d^2) * t(svdXtilde$v))))
      s <- sdBetas/min(sdBetas) #covariate specific posterior sd divided by minimum sd
      
      output<-lapply(maxsel2,function(x){
        fsel <- function(t,maxselec = x){
          An <- which(abs(beta[ind]) > t*s) #marginal posterior credible set selection rule for some t
          
          return(length(An)- maxselec) #number of non-zero penalised covariates
        }
        maxT <- max(abs(beta[ind])/s)*2
        t <- uniroot(fsel, interval = c(0, maxT), maxiter = 200,tol=10^(-20))$root
        indPost <- ind[which(abs(beta[ind]) > t*s)] #index of selected penalised covariates
        indPostpen <- which(which(penfctr!=0) %in% indPost) #index in 1:length(pen) of selected penalised covariates
        indAll <- c(indPost,which(penfctr==0)) #index of selected penalised covariates and unpenalised covariates

        
        if(grepl("dense",postselection)){
          if(x==maxsel2[1]) if(!silent) print("Selected covariates are refit with previously estimated weighted ridge prior")
          
          lambdaoverall <- sigmahat/tautrgt
          lam2 <- sigmahat/tautrgt/n*sd_y
          #Recompute beta using only selected beta and group ridge penalty (without lasso penalty)
          offset <- rep(0,n)
          if(length(muhatp)==1) muhatp <- rep(muhatp,sum(penfctr!=0))
          if(!all(muhatp==0)) offset <- X[,indAll, drop=F] %*% muhatp[indPostpen, drop=F]
          Xacc <- X
          Xacc[,pen] <- as.matrix(X[,pen] %*% sparseMatrix(i=1:length(lambdap),j=1:length(lambdap),
                                                           x=c(1/sqrt(lambdap/lambdaoverall))) )

          if(model=="cox"){
            glmPost <- glmnet(Xacc[,indAll, drop=F],Y,alpha=0,
                            lambda = lam2,family=fml,
                            offset = offset ,standardize = F,
                            penalty.factor=penfctr[indAll], thresh=10^-10)
          }else{
            glmPost <- glmnet(Xacc[,indAll, drop=F],Y,alpha=0,
                            lambda = lam2,family=fml,
                            offset = offset , intercept = intrcpt, standardize = F,
                            penalty.factor=penfctr[indAll], thresh=10^-10)
          }
          betaPost <- rep(0,p)
          betaPost[indAll] <- as.vector(glmPost$beta)
          betaPost[indPost] <- c(1/sqrt(lambdap[indPostpen]/lambdaoverall)) * betaPost[indPost] + muhatp[indPostpen]
        }else{ #sparse; refit with newly cross-validated lambda
          if(model=="cox"){
            lambdaGLM<-cv.glmnet(X[,indAll, drop=F],Y,alpha=0,family=fml,
                                 standardize = F,penalty.factor=penfctr[indAll]) #alpha=0 for ridge
            lam2<-lambdaGLM$lambda.min
          }else{
            lambdaGLM<-cv.glmnet(X[,indAll, drop=F],Y,alpha=0,family=fml,
                                 standardize = F,penalty.factor=penfctr[indAll]) #alpha=0 for ridge
            lam2<-lambdaGLM$lambda.min
          }
          if(model=="cox"){
            glmPost <- glmnet(X[,indAll, drop=F],Y,alpha=0,
                            lambda = lam2,family=fml,
                            offset = offset ,standardize = F,
                            penalty.factor=penfctr[indAll], thresh=10^-10)
          }else{
            glmPost <- glmnet(X[,indAll, drop=F],Y,alpha=0,
                            lambda = lam2,family=fml,
                            intercept = intrcpt, standardize = F,
                            penalty.factor=penfctr[indAll], thresh=10^-10)
          }
          betaPost <- rep(0,p)
          betaPost[indAll] <- as.vector(glmPost$beta)
        }
        
        whichPost <- indPost
        output<-list()
        output$betaPost <- betaPost #post-selection beta using Elastic Net
        output$whichPost <- whichPost #index of selected covariates
        output$a0 <- glmPost$a0
        #output$offsetPost <- offset #offset used in Post

        if(!all(is.nan(X2))){
          if(model=="linear"){
            YpredPost <- X2 %*% c(betaPost) + glmPost$a0
            MSEPost <- sum((YpredPost-Y2)^2)/length(Y2)
          } 
          if(model=='logistic'){
            YpredPost <- 1/(1+exp(-(X2 %*% c(betaPost) + glmPost$a0)))
            MSEPost<- sum((YpredPost-Y2)^2)/length(Y2)
            if(any(is.nan(YpredPost))){browser()}
          }else if(model=='cox'){ 
            expXbPost<-exp(X %*% c(betaPost))
            h0 <- sapply(1:length(Y[,1]),function(i){Y[i,2]/sum(expXbPost[Y[,1]>=Y[i,1]])})#updated baseline hazard in censored times for left out samples
            H0 <- sapply(Y2[,1],function(Ti){sum(h0[Y[,1]<=Ti])})
            YpredPost <- H0*exp(X2 %*% c(betaPost))
            MSEPost<- sum((YpredPost-Y2[,2])^2)/length(Y2[,2])
          }
          
          output$MSEPost <- MSEPost #MSE on independent data set (if given)
          output$YpredPost <- YpredPost #predictions for independent data set (if given)
        }
        
        return(output)
      })
      
    }
    
  }
  
  #reshape output data
  res<-list()
  size <- length(maxsel2)
  for(attr in names(output[[1]])){
    if(attr != "whichPost"){
      res[[attr]] <- matrix(unlist(lapply(output,function(x){x[[attr]]})),length(output[[1]][[attr]]),size)
      colnames(res[[attr]])<-paste("MaxSelec: ",maxsel2,sep="")
    }  
  }
  return(res)
}

#Produce balanced folds----
produceFolds <- function(nsam,outerfold,response,model="logistic",balance=TRUE,fixedfolds=F){
  if(fixedfolds) set.seed(3648310) #else set.seed(NULL)
  if(model=="linear") balance=F
  if(!balance){
    rand<-sample(1:nsam)
    grs1 <- floor(nsam/outerfold)
    grs2 <- grs1+1
    ngr1 <- outerfold*grs2 - nsam
    folds <- lapply(1:outerfold,function(xg) {
      if(xg <= ngr1) els <- rand[(1+(xg-1)*grs1):(xg*grs1)] else els <- rand[(ngr1*grs1 + 1+(xg-ngr1-1)*grs2):(ngr1*grs1 + (xg-ngr1)*grs2)]
      return(els)
    }
    )} else {  #balanced folds
      if(model=="logistic") if(class(response)=="factor") nev <- which((as.numeric(response)-1)==1) else nev <- which(response==1)  
      if(model=="survival") nev <- which(response[,1]==1)    
      nsamev <- length(nev) 
      randev<-sample(nev)
      grs1 <- floor(nsamev/outerfold)
      grs2 <- grs1+1
      ngr1 <- outerfold*grs2 - nsamev
      foldsev <- lapply(1:outerfold,function(xg) {
        if(xg <= ngr1) els <- randev[(1+(xg-1)*grs1):(xg*grs1)] else els <- randev[(ngr1*grs1 + 1+(xg-ngr1-1)*grs2):(ngr1*grs1 + (xg-ngr1)*grs2)]
        return(els)
      }
      )
      nonev <- setdiff(1:nsam,nev)
      nsamnonev <- length(nonev) 
      randnonev<-sample(nonev)
      grs1 <- floor(nsamnonev/outerfold)
      grs2 <- grs1+1
      ngr1 <- outerfold*grs2 - nsamnonev
      foldsnonev <- lapply(1:outerfold,function(xg) {
        if(xg <= ngr1) els <- randnonev[(1+(xg-1)*grs1):(xg*grs1)] else els <- randnonev[(ngr1*grs1 + 1+(xg-ngr1-1)*grs2):(ngr1*grs1 + (xg-ngr1)*grs2)]
        return(els)
      }
      )
      folds <- lapply(1:outerfold,function(i) c(foldsev[[i]],foldsnonev[[i]]))
    }
  return(folds)
}

#Estimate maximum marginal likelihood estimates for linear model----
mlestlin <- function(Y,X,lambda=NaN,sigmasq=NaN,mu=NaN,tausq=NaN){
  #lambda,sigmasq,mu are possibly fixed
  maxv <- var(Y)
  p<-dim(X)[2]
  n<-dim(X)[1]
  prms <- paste(c("l","s","m","t")[is.nan(c(lambda,sigmasq,mu,tausq))],collapse='') #unknown parameters: l lambda, s sigma, m mu
  if(prms=='t'||prms=='mt'){tausq <- sigmasq/lambda; prms <- gsub("t","",prms)}
  if(prms=='l'||prms=='lm'){lambda <- sigmasq/tausq; prms <- gsub("l","",prms)}
  if(prms=='s'||prms=='sm'){sigmasq <- lambda*tausq; prms <- gsub("s","",prms)}
  switch(prms,
         'lsmt'={ #lambda, sigma, mu, tau unknown
           sim2 = function(ts){
             tausq<-ts[1];sigmasq<-ts[2];mu<-ts[3]
             n <- nrow(X)
             varY <- X %*% t(X) * exp(tausq) + diag(rep(1,n))*exp(sigmasq)
             meanY <- X%*%rep(mu,p)
             mlk <- -dmvnorm(c(Y),mean=meanY,sigma=varY,log=T)
             return(mlk)
           }
           op <- optim(c(log(0.01),log(maxv),0),sim2)
           tausq <- exp(op$par[1]); sigmasq <- exp(op$par[2])
           mu <- op$par[3]; lambda <- sigmasq/tausq
         },
         'lsm'={ #lambda, sigma, mu unknown, tau known
           sim2 = function(ts){
             sigmasq<-ts[1];mu<-ts[2]
             n <- nrow(X)
             varY <- X %*% t(X) * exp(tausq) + diag(rep(1,n))*exp(sigmasq)
             meanY <- X%*%rep(mu,p)
             mlk <- -dmvnorm(c(Y),mean=meanY,sigma=varY,log=T)
             return(mlk)
           }
           op <- optim(c(log(maxv),0),sim2)
           sigmasq <- exp(op$par[1])
           mu <- op$par[2]; lambda <- sigmasq/tausq
         },
         'lst'={ #lambda, sigma, tau unknown, mu known
           sim2 = function(ts){
             tausq<-ts[1];sigmasq<-ts[2]
             n <- nrow(X)
             varY <- X %*% t(X) * exp(tausq) + diag(rep(1,n))*exp(sigmasq)
             meanY <- X%*%rep(mu,p)
             mlk <- -dmvnorm(c(Y),mean=meanY,sigma=varY,log=T)
             return(mlk)
           }
           op <- optim(c(log(0.01),log(maxv)),sim2)
           tausq <- exp(op$par[1]); sigmasq <- exp(op$par[2])
           lambda <- sigmasq/tausq
         },
         'lmt'={ #lambda, tau, mu unknown, sigma known
           sim2 = function(ts){
             tausq<-ts[1]; mu <- ts[2]
             n <- nrow(X)
             varY <- X %*% t(X) * exp(tausq) + diag(rep(1,n))*sigmasq
             meanY <- X%*%rep(mu,p)
             mlk <- -dmvnorm(c(Y),mean=meanY,sigma=varY,log=T)
             return(mlk)
           }
           op <- optim(c(log(0.01),0),sim2)
           tausq <- exp(op$par[1])
           mu <- op$par[2]; lambda <- sigmasq/tausq
         },
         'lt'={ #lambda, tau unknown, sigma, mu known
           sim2 = function(ts){
             tausq<-ts[1];
             n <- nrow(X)
             varY <- X %*% t(X) * exp(tausq) + diag(rep(1,n))*sigmasq
             meanY <- X%*%rep(mu,p)
             mlk <- -dmvnorm(c(Y),mean=meanY,sigma=varY,log=T)
             return(mlk)
           }
           op <- optim(c(log(0.01)),sim2,method="Brent",lower=log(1e-10),upper=log(100))
           #op <- optimize(sim2,c(-100,100))
           #op <- optim(c(log(0.01)),sim2)
           #browser()
           tausq <- exp(op$par[1])
           lambda <- sigmasq/tausq
         },
         'smt'={ #sigma, tau, mu unknown, lambda known
           sim2 = function(ts){
             tausq<-log(exp(ts[1])/lambda); sigmasq<-ts[1]; mu<-ts[2]
             n<- nrow(X)
             varY <- X %*% t(X) * exp(tausq) + diag(rep(1,n))*exp(sigmasq)
             meanY <- X%*%rep(mu,p)
             mlk <- -dmvnorm(c(Y),mean=meanY,sigma=varY,log=T)
             return(mlk)
           }
           op <- optim(c(log(maxv),0),sim2,method="L-BFGS-B",lower=c(log(1e-10),-10^{99}),upper=c(log(2*maxv),10^{99}))
           tausq <- exp(op$par[1])/lambda; sigmasq <- exp(op$par[1])
           mu <- op$par[2]
           #or for small n
           #sigmasq <- sum(diag(solve(X%*%t(X)+lambda*diag(n),lambda*Y%*%t(Y))))/n
           #mlests <- c(sigmasq,sigmasq/lambda)
         },
         'st'={ #sigma, tau unknown, lambda, mu known
           sim2 = function(ts){
             tausq<-log(exp(ts)/lambda); sigmasq<-ts
             n<- nrow(X)
             varY <- X %*% t(X) * exp(tausq) + diag(rep(1,n))*exp(sigmasq)
             meanY <- X%*%rep(mu,p)
             mlk <- -dmvnorm(c(Y),mean=meanY,sigma=varY,log=T)
             return(mlk)
           }
           op <- optim(c(log(maxv)),sim2,method="Brent",lower=log(1e-10),upper=log(2*maxv))
           tausq <- exp(op$par[1])/lambda; sigmasq <- exp(op$par[1])
         },
         'ls'={ #lambda, sigma unknown, tau, mu known
           sim2 = function(ts){
             sigmasq<-ts[1];
             n <- nrow(X)
             varY <- X %*% t(X) * exp(tausq) + diag(rep(1,n))*sigmasq
             meanY <- X%*%rep(mu,p)
             mlk <- -dmvnorm(c(Y),mean=meanY,sigma=varY,log=T)
             return(mlk)
           }
           op <- optim(c(log(maxv)),sim2,method="Brent",lower=log(1e-10),upper=log(2*maxv))
           sigmasq <- exp(op$par[1])
           lambda <- sigmasq/tausq
         },
         'm'={ #mean unknown, lambda, sigma, tau known
           sim2 = function(ts){
             mu<-ts
             n<- nrow(X)
             varY <- X %*% t(X) * tausq + diag(rep(1,n))*sigmasq
             meanY <- X%*%rep(mu,p)
             mlk <- -dmvnorm(c(Y),mean=meanY,sigma=varY,log=T)
             return(mlk)
           }
           op <- optim(0,sim2,method="Brent",lower=-10^{99},upper=10^{99})
           mu <- op$par[1]
         }
  )
  mlests <- c(lambda,sigmasq,mu,tausq)
  return(mlests)
}



#simulate data for linear or logistic----
simDat <- function(n,p,n2=20,muGrp,varGrp,indT,sigma=1,model='linear',flag=F){
  #Simulate data:
  #X~MVN(0,1)
  #beta_k~N(mu_{g(k)},sig_{g(k)}
  #Linear regression: Y~N(X*beta,sigma^2)
  #Logistic regression: Y~binom(expit(X*beta))
  #
  #Input:
  #n: number of observations
  #p: number of covariates
  #n2: number of observations of independent data sample
  #muGrp: prior mean of grouped betas
  #varGrp: variance of grouped betas
  #indT: px1 vector of group indices 1,..,G
  #model: generate data for "linear" or "logistic" setting
  #flag: True if plots and histograms of generated data have to be shown
  
  #generate expression data
  X <- rmvnorm(n,rep(0,p),diag(rep(1,p))) #X~MVN(0,1)
  X2 <- rmvnorm(n2,rep(0,p),diag(rep(1,p))) #X2 independent data
  
  #center X
  meansX<-apply(X,2,mean)
  Xctd <- t(t(X) - meansX) #center X
  meansX2<-apply(X2,2,mean)
  X2ctd <- t(t(X2) - meansX2) #center X
  beta <- c(rmvnorm(1,muGrp[indT],
                  diag(varGrp[indT]))) #beta_k~N(mu_{g(k)},tau_{g(k)}^2)
  
  if(model=='linear'){
    Y <- Xctd %*% beta + rnorm(n,0,sd=sigma) #Y~N(X*beta,sigma^2)
    Y2 <- X2ctd %*% beta + rnorm(n2,0,sd=sigma) #Y~N(X*beta,sigma^2)
    
    if(flag){
      if(!silent) print(paste("Simulated data for",model,"regression"))
      #plot data
      plot(Y)
      points(Xctd %*% beta, col='red') #Y without the added noise
      legend('right',c('Y','X*beta'),pch=c(1,1),col=c('black','red'))
      hist(Y)
    }
    
  } 
  if(model=='logistic'){
    expXb <- exp(Xctd %*% beta)
    Y <-  rbinom(n,1,expXb/(1+expXb)) #Y~binom(expit(X*beta))
    expX2b <- exp(X2ctd %*% beta)
    Y2 <- rbinom(n2,1,expX2b/(1+expX2b))
    
    if(flag){
      if(!silent) print(paste("Simulated data for",model,"regression"))
      #plot data
      plot(Y)
      points(expXb/(1+expXb), col='red') #Y without the added noise
      legend('right',c('Y','expit(X*beta)'),pch=c(1,1),col=c('black','red'))
      hist(expXb/(1+expXb))
    }
  }
  
  return(
    list(
      beta=beta,
      Xctd=Xctd, 
      Y=Y,
      X2ctd=X2ctd,
      Y2=Y2
    )
  )
}

#Perform cross-validation----
cv.ecpc <- function(Y,X,type.measure="MSE",outerfolds=10,
                    lambdas=NULL,ncores=1,balance=T,silent=F,...){
  ecpc.args <- list(...)
  if(!is.element("model",names(ecpc.args))){
    if(all(is.element(Y,c(0,1))) || is.factor(Y)){
      model <- "logistic" 
    } else if(all(is.numeric(Y)) & !(is.matrix(Y) && dim(Y)[2]==2)){
      model <- "linear"
    }else{
      model <- "cox"
    }
  }else{
    model <- ecpc.args$model; ecpc.args$model <- NULL
  }
  if(!is.element("postselection",names(ecpc.args))){
    postselection <- "elnet,dense"
  }else{
    postselection <- ecpc.args$postselection; ecpc.args$postselection <- NULL
  }
  if(!is.element("maxsel",names(ecpc.args))){
    maxsel <- 10
  }else{
    maxsel <- ecpc.args$maxsel; ecpc.args$maxsel <- NULL
  }
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  if(is.numeric(outerfolds)){
    folds2<-produceFolds(n,outerfolds,Y,balance=balance,model=model) #produce folds balanced in response
  }else{
    folds2 <- outerfolds
  }
  nfolds <- length(folds2)
  if(length(lambdas)==0){
    lambdas <- rep(NaN,length(folds2))
  }
  
  grpsno <- c(unlist(sapply(ecpc.args$groupings,function(x){1:length(x)}))) #vector with group numbers in all groupings
  grpngsno <- c(unlist(sapply(1:length(ecpc.args$groupings),function(i){rep(i,length(ecpc.args$groupings[[i]]))}))) #vector with groupings numbers
  dfGrps<-data.frame() #data frame in which group and grouping weights are stored
  df<-data.frame() #data frame in which predicted values are stores for each of the samples in the left-out fold
  Res<-list() #list in which raw output of ecpc is stored (e.g. estimated regression coefficients)
  if(ncores>1){
    cl <- makeCluster(ncores) #set up parallel cluster
    registerDoParallel(cl)
    finalMatrix <- foreach(i=1:nfolds, .combine=rbind, 
                           .packages = c("glmnet","penalized","mvtnorm","gglasso",
                                         "Matrix","Rsolnp","ecpc")) %dopar% {
         
         tic<-proc.time()[[3]]
         Res[[i]]<-do.call(ecpc,args=c(list(Y=Y[-folds2[[i]]],X=X[-folds2[[i]],],Y2=Y[folds2[[i]]],
                                            X2=X[folds2[[i]],],lambda=lambdas[i],postselection=postselection,
                                            maxsel = maxsel,model=model,silent=silent),ecpc.args))
         Res[[i]]$time <- proc.time()[[3]]-tic
         
         if(postselection!=F){
           df2<-data.frame("Ypred"=c(Res[[i]]$Ypred,Res[[i]]$Ypredridge,c(Res[[i]]$YpredPost)))
           df2$Method <- rep(c("ecpc","ordinary.ridge",paste("ecpc",maxsel,"vars",sep="")),each=length(folds2[[i]]))
           df2$NumberSelectedVars <- rep(c(sum(Res[[i]]$beta!=0),p,maxsel),each=length(folds2[[i]]))
           df2$Fold <- i
           df2$Sample <- rep(folds2[[i]],2+length(maxsel))
           df2$Time <-  Res[[i]]$time
           df2$Truth <- rep(Y[folds2[[i]]],2+length(maxsel))
         }else{
           df2<-data.frame("Ypred"=c(Res[[i]]$Ypred,Res[[i]]$Ypredridge))
           df2$Method <- rep(c("ecpc","ordinary.ridge"),each=length(folds2[[i]]))
           df2$NumberSelectedVars <- rep(c(sum(Res[[i]]$beta!=0),p),each=length(folds2[[i]]))
           df2$Fold <- i
           df2$Sample <- rep(folds2[[i]],2)
           df2$Time <-  Res[[i]]$time
           df2$Truth <- rep(Y[folds2[[i]]],2)
         }
         
         df3<-data.frame("Group"=grpsno,
                         "Grouping"=grpngsno,
                         "Group weight"=Res[[i]]$gamma,
                         "Grouping weight"=Res[[i]]$w[grpngsno])
         df3$Tau.ecpc <- Res[[i]]$tauglobal #global tau^2
         df3$Tau.ridge <- 1/Res[[i]]$lambdaridge #ordinary ridge tau^2
         df3$Method <- "ecpc"
         df3$Fold <- i
         
         if(!silent) print(paste(Sys.time(),"fold",i,"of",nfolds,"done"))
         
         list("Res"=Res,"df"=df2,"dfGrps"=df3)
       }
    
    Res <- lapply(1:nfolds,function(i) finalMatrix[i,1][[1]][[i]])
    df2 <- lapply(1:nfolds,function(i) finalMatrix[i,2][[1]])
    dfGrps2 <- lapply(1:nfolds,function(i) finalMatrix[i,3][[1]])
    df <- df2[[1]]; for(i in 2:nfolds) df <- rbind(df,df2[[i]])
    dfGrps <- dfGrps2[[1]]; for(i in 2:nfolds) dfGrps <- rbind(dfGrps,dfGrps2[[i]])
    stopCluster(cl); rm(cl)
  }else{
    for(i in 1:nfolds){
       tic<-proc.time()[[3]]
       Res[[i]]<-do.call(ecpc,args=c(list(Y=Y[-folds2[[i]]],X=X[-folds2[[i]],],Y2=Y[folds2[[i]]],
                                          X2=X[folds2[[i]],],lambda=lambdas[i],postselection=postselection,
                                          maxsel = maxsel,model=model,silent=silent),ecpc.args))
       Res[[i]]$time <- proc.time()[[3]]-tic
       if(postselection!=F){
         df2<-data.frame("Ypred"=c(Res[[i]]$Ypred,Res[[i]]$Ypredridge,c(Res[[i]]$YpredPost)))
         df2$Method <- rep(c("ecpc","ordinary.ridge",paste("ecpc",maxsel,"vars",sep="")),each=length(folds2[[i]]))
         df2$NumberSelectedVars <- rep(c(sum(Res[[i]]$beta!=0),p,maxsel),each=length(folds2[[i]]))
         df2$Fold <- i
         df2$Sample <- rep(folds2[[i]],2+length(maxsel))
         df2$Time <-  Res[[i]]$time
         df2$Truth <- rep(Y[folds2[[i]]],2+length(maxsel))
       }else{
         df2<-data.frame("Ypred"=c(Res[[i]]$Ypred,Res[[i]]$Ypredridge))
         df2$Method <- rep(c("ecpc","ordinary.ridge"),each=length(folds2[[i]]))
         df2$NumberSelectedVars <- rep(c(sum(Res[[i]]$beta!=0),p),each=length(folds2[[i]]))
         df2$Fold <- i
         df2$Sample <- rep(folds2[[i]],2)
         df2$Time <-  Res[[i]]$time
         df2$Truth <- rep(Y[folds2[[i]]],2)
       }
       df<-rbind(df,df2)
       
       df3<-data.frame("Group"=grpsno,
                       "Grouping"=grpngsno,
                       "Group weight"=Res[[i]]$gamma,
                       "Grouping weight"=Res[[i]]$w[grpngsno])
       df3$Tau.ecpc <- Res[[i]]$tauglobal #global tau^2
       df3$Tau.ridge <- 1/Res[[i]]$lambdaridge #ordinary ridge tau^2
       df3$Method <- "ecpc"
       df3$Fold <- i
       dfGrps<-rbind(dfGrps,df3)
       
       if(!silent) print(paste(Sys.time(),"fold",i,"of",nfolds,"done"))
    }
  }

  df$Method<-as.factor(df$Method)
  dfGrps$Group<-as.factor(dfGrps$Group)
  
  #data frame with performance measure
  if(is.factor(df$Truth)){
    warning("Response Y given as factor, transformed to numeric to compute AUC")
    if(!silent) print(levels(df$Truth)[1],"transformed to",0)
    if(!silent) print(levels(df$Truth)[2],"transformed to",1)
    df$Truth <- as.numeric(df$Truth)-1
  }
  if(type.measure=="MSE"){
    dfCVM <- df %>% group_by(Method,Fold) %>% summarise(CVM = mean((Ypred-Truth)^2),Type="MSE",
                                                   NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  }
  else if(type.measure=="AUC"){
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
    dfCVM <- dfROC %>% group_by(Method) %>% summarise(CVM=mean(AUC),Type="AUC",
                                                      NumberSelectedVars=mean(NumberSelectedVars)) %>% ungroup()
  }else{
    warning(paste("The type of measure",type.measure,"is not yet supported."))
  }
  
  return(list(ecpc.fit=Res, #list with model fit in each fold
              dfPred = df, #data frame with information about out-of-bag predictions
              dfGrps=dfGrps, #data frame with information about estimated group and grouping weights across folds
              dfCVM = dfCVM #data frame with cross-validated performance metric
  ))
}

#Discretise continuous co-data hierarchically----
#Output: list of groups of covariates of varying size
splitMedian <- function(values,index,depth,minGroupSize=50,first=T,split="both"){
  #split: "both" or "lower" if both groups have to be split recursively or only lower, lower-valued group
  if(missing(depth)){
    depth <- floor(log2(length(values)/minGroupSize)) #number of layers in hierarchical tree such that lowest level contains at least minSize group members (start from 2 groups)
  }
  if(missing(index)&first){
    index <- 1:length(values)
  }
  
  medianX <- median(values) #compute median
  indSmall <- values<=medianX #index of values smaller than median
  indLarge <- values>medianX #index of values larger than median

  if(first){
    groups <- list(index,index[indSmall],index[indLarge]) #split group at median, include group with all variables
  }else{
    groups <- list(index[indSmall],index[indLarge]) #split group at median
  }
  if(depth>1){ #TD: improve recursion for leaves only splitting in one side
    if(split=="both"){
      return(c(groups,
               splitMedian(values[indSmall],index[indSmall],depth-1,minGroupSize=minGroupSize,first=F,split=split),
               splitMedian(values[indLarge],index[indLarge],depth-1,minGroupSize=minGroupSize,first=F,split=split)))
    }else if(split=="lower"){
      return(c(groups,
               splitMedian(values[indSmall],index[indSmall],depth-1,minGroupSize=minGroupSize,first=F,split=split)))
    }else{
      stop("split should be one of both or lower")
    }
  }else{
    #Stop if one of the two groups has less than minGroupSize variables (could happen if lots of variables have the same value)
    LargeEnough <- sapply(groups,function(x){length(x)>=minGroupSize})
    if(all(LargeEnough)){
      return(groups)
    }else{
      return(NULL)
    }
  }
}

#Obtain grouping on group level for hierarhical groups----
#Output: list of groups of group numbers (i.e. on group level) defining the hierarchy
obtainHierarchy <- function(grouping,penalty="LOG"){
  if(penalty=="LOG"){
    #if Latent Overlapping Group Lasso hypershrinkage is used (for now the only option in ecpc), 
    # the hierarchy on the covariate groups is defined as:
    # for each group number (node in the hierarchical tree),
    # make a group of the group number, say g, and all group numbers of groups that are supersets of group g
    grouping.grplvl <- lapply(grouping,function(i){as.numeric(which(sapply(grouping,function(j){all(i%in%j)})))})
  }else if(penalty=="GL"){
    #NOTE: this option is not yet available in ecpc, use LOG for hierarchical groups;
    #if Group Lasso hypershrinkage is used, the hierarchy on the covariate groups is defined as:
    # for each group number (node in the hierarchical tree): 
    # make a group of the group number, say g, and all group numbers of groups that are subsets of group g
    grouping.grplvl <- lapply(grouping,function(j){as.numeric(which(sapply(grouping,function(i){all(i%in%j)})))})
  }
  return(grouping.grplvl)
}

#Fit hierarchical lasso using LOG penalty----
hierarchicalLasso <- function(X,Y,grouping){
  #X: nxp matrix with observed data
  #Y: nx1 vector with response
  #grouping: list with hierarchical group indices
  
  p<-dim(X)[2]
  Xxtnd <- X[,unlist(grouping)] #extend matrix such to create artifical non-overlapping groups
  #create new group indices for Xxtnd
  Kg2 <- c(1,sapply(grouping,length)) #group sizes on group level (1 added to easily compute hier. group numbers)
  G2 <- length(Kg2)-1
  groupxtnd <- lapply(2:length(Kg2),function(i){sum(Kg2[1:(i-1)]):(sum(Kg2[1:i])-1)}) #list of indices in each group
  groupxtnd2 <- unlist(sapply(1:G2,function(x){rep(x,Kg2[x+1])})) #vector with group number
  
  cvfit <- cv.gglasso(x=Xxtnd,y=Y, group = groupxtnd2,pred.loss = c("L2"))
  
  #Hierarchical group lasso estimate for group variances
  #fit2<-gglasso(x=Xxtnd,y=Y,group = groupxtnd2, loss="ls")
  intrcptOG <- coef(cvfit,s=cvfit$lambda.min)[1]
  vtilde <- coef(cvfit,s=cvfit$lambda.min)[-1]
  v<-lapply(groupxtnd,function(g){
    x<-rep(0,p)
    x[unlist(grouping)[g]]<-x[unlist(grouping)[g]]+vtilde[g]
    return(x)
  })
  betahat <- c(apply(array(unlist(v),c(p,G2)),1,sum))
  return(list("betas"=betahat,"a0"=intrcptOG, 
              "lambda"=cvfit$lambda.min,"group.weights"=sqrt(Kg2[-1])))
}

#Multiply slices of 3D-arrays----
.matmult3d <- function(A,B){
  dimA <- dim(A)
  dimB <- dim(B)
  if(length(dimB)==2){                      #B matrix not an array, A 3D array
    return(array(apply(A,3,function(x){x%*%B}),c(dimA[1],dimB[2],dimA[3])))
  } else if(length(dimA)==2){               #A matrix not an array, B 3D array
    return(array(apply(B,3,function(x){A%*%x}),c(dimA[1],dimB[2],dimB[3])))
  } else if(dimB[3]==dimA[3]){              #A*B "matrix-wise" nsplits different A with nsplits different B
    return(array(apply(array(1:dimB[3],c(dimB[3],1,1)),1,function(i){A[,,i]%*%B[,,i]}),c(dimA[1],dimB[2],dimB[3])))
  } else if(dimB[3]==1){                    #A*B for nsplits different A and one B
    return(array(apply(A,3,function(x){x%*%B[,,1]}),c(dimA[1],dimB[2],dimA[3])))
  } else if(dimA[3]==1){                    #A*B for one A and nsplits different B
    return(array(apply(B,3,function(x){A[,,1]%*%x}),c(dimA[1],dimB[2],dimB[3])))
  } else{warning("matmult3d error")
    browser()}
}

#Iterative weighted least squares algorithm for logistic regression----
.IWLS <- function(X,Y,penalties,targets,eps=1e-7,maxItr=100){
  #Input:
  #X: nxp data matrix
  #Y: nx1 response
  #penalties: px1 vector with penalties (prop. to prior variance)
  #targets (optional): px1 vector with targets (equiv. to prior mean)
  #eps: numerical bound for convergence
  #maxItr: maximum number of iterations used in IWLS
  #Output:
  #betas: estimate for regression coefficients
  #convergence (T/F): T if IWLS has converged under threshold eps, F otherwise 
  #nIt: number of iterations needed for convergence
  
  if(missing(targets)) targets <- rep(0,dim(X)[2])
  
  betaold <- rep(0,p) #intial value
  it <- 0; nrm <- Inf
  while(nrm>eps & it<maxItr){
    Xb <- LogData$Xstd %*% betaold #linear predictor
    Ypred<-1/(1+exp(-Xb)) #predicted probability 
    W<-diag(c((Ypred*(1-Ypred)))) #Var(Y)
    
    XtXDinv <- solve(t(LogData$Xstd) %*% W %*% LogData$Xstd + 2*diag(lambda,p)) #inverting pxp matrix
    z <-  Xb + (LogData$Y - Ypred)/c(diag(W))
    betas <- XtXDinv %*% (t(LogData$Xstd) %*% W %*% z + 2*diag(lambda,p)%*%mup)
    nrm <- sum((betas-betaold)^2)
    betaold <- betas
    it <- it + 1
  }
  convergence <- nrm<eps
  return(list(betas=betas,convergence=convergence,nIt=it))
}

#Compute inverse gamma penalty parameters with mode 1 for given group size and hyperpenalty lambda----
prmsIGMode1 <- function(lam,Kg){
  #Input:
  #lam: penalty parameter lambda\in(0,\infty)
  #Kg: vector of group sizes
  #Output:
  #list with vector of alpIG, betIG:
  #the parameters of inverse gamma penalty such that:
  #Mode=1: betIG/(alpIG+1)=1
  #Variance=1/(lam*Kg): betIG^2/(alp-1)^2/(alp-2)=1/(lam*Kg)
  
  const <- 1/(lam*Kg)
  alpIG <- sapply(const,function(const){
    #solve cubic equation
    a <- 1
    b <- -(4+1/const)
    c <- 5-2/const
    d<- -(2+1/const)
    
    Delt0 <- b^2 - 3*a*c
    Delt1 <- 2*b^3 - 9*a*b*c + 27*a^2*d
    C <- (0.5*(Delt1+sqrt(as.complex(Delt1^2 - 4*Delt0^3))))^(1/3)
    if(Re(C)==0&&Im(C)==0) C <- (0.5*(Delt1-sqrt(Delt1^2 - 4*Delt0^3)))^(1/3)
    
    xi <- (0.5*(-1+sqrt(as.complex(-3))))^c(0,1,2)
    
    x <- -1/(3*a)*(b + xi*C + Delt0/(xi * C)) 
    xReal <- Re(x[abs(Im(x))<10^-12]) #real roots
    if(length(xReal)==0) xReal <- Re(x[which.min(abs(Im(x)))])
    return(xReal)
  })
  betIG <- alpIG+1
  
  return(list(alpIG,betIG))
}

###Plotting functions
#Visualise grouping----
visualiseGrouping <- function(Grouping,groupweights,grouping.grouplvl,nodeSize=10,ls=1){
  #return ggplot object with graphical visualisation of one grouping:
  #graph with nodes, possibly connected if hierarchy is given
  #if group weights are given, nodes are coloured accordingly: 
  #  (blue if >1, red if <1, ivory if =1 or no group weight given, white if not selected (=0))
  #
  #Input:
  #Grouping: list of G elements containing the indices of the variables in that group
  #grouping.grouplvl (optional): hierarchical structure on the groups
  #groupweights (optional): group weights
  #Output:
  #Plot in a ggplot2 object
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cols<-c("#FC8D59","ivory","#91BFDB") #brewer.pal(3,"RdYlBu")[c(1,3)]
  
  Kg <- sapply(Grouping,length) #group size of each group
  G<-length(Grouping) #number of groups in grouping
  groupNames <- names(Grouping) #group names
  if(length(names(Grouping))==0) groupNames <- 1:G
  if(missing(grouping.grouplvl)||length(grouping.grouplvl)==0){
    basicGraph <- T
    grouping.grouplvl <- lapply(1:G,function(x) x)
  }else basicGraph <- F
  i<-unlist(sapply(1:G,function(x){rep(x,unlist(Kg)[x])}))
  j<-unlist(unlist(Grouping))
  ind <- sparseMatrix(i,j,x=1) #sparse matrix with ij element 1 if jth element in group i (global index), otherwise 0
  
  #Co-data matrix: sparse matrix with ij element 1/Ij if beta_j in group i
  Zt<-ind; 
  if(G[1]>1){
    Zt[1:G[1],]<-t(t(ind[1:G[1],])/apply(ind[1:G[1],],2,sum))
  }
  
  #create data frame with nodes
  # We can add a second data frame with information for each node!
  vertices <- data.frame(name=groupNames)
  if(missing(groupweights)){
    vertices <- vertices %>% mutate("name2"=name,"Selected"=T,"GrpWeight"=rep(1,G))
    showFillLegend <- F
  } 
  else{
    vertices <- vertices %>% mutate("GrpWeight"=round(groupweights,digits=2),
                                    "Selected"=groupweights!=0,
                                    "name2"=paste(name," (",GrpWeight,")",sep=""))
    showFillLegend <- T
  }  
  
  temp<-obtainHierarchy(grouping.grouplvl)
  children<-lapply(1:length(temp),function(i){
    children <- setdiff(temp[[i]],i)
    for(j in children){
      children <- setdiff(children,setdiff(temp[[j]],j))
    }
    if(length(children)==0){return(NULL)}
    return(children)
  })
  # create an edge list data frame giving the hierarchical structure of your individuals
  edgeList <- data.frame(from=as.character(unlist(sapply(1:length(children),function(x){rep(x,length(children[[x]]))}))),
                         to=as.character(unlist(children)))
  mygraph <- graph_from_data_frame( edgeList ,vertices=vertices)
  
  if(basicGraph){ # Basic groups (no hierarchy)
    # Create a graph object
    mygraph <- graph_from_data_frame( edgeList ,vertices=vertices)
    
    lay <- create_layout(mygraph, layout = 'linear')
    p<-ggraph(lay) +
      geom_node_label(aes(label=name2,fill=GrpWeight,col=Selected), show.legend = showFillLegend, fontface = "bold",
                      size=nodeSize,repel=T,
                      point.padding = NA, box.padding = 0, force = 0.1)+
      #scale_fill_gradientn(colors=c("white",cols),breaks=c(0,1,NA),limits=c(0,NA))+
      scale_fill_gradientn(colors=c("white",cols),
                           #breaks=c(0,1,NA),
                           limits=c(0,NA),
                           values = scales::rescale(
                             c(0,min(10^-6,lay$GrpWeight[lay$GrpWeight>0]), #white for 0
                               min(10^-6,lay$GrpWeight[lay$GrpWeight>0]),1-1e-6,  #red for <1
                               1-1e-6,1+1e-6,
                               1+1e-6, max(lay$GrpWeight))),
                           name="Group weight")+
      scale_color_manual(values=c("TRUE"="black","FALSE"=cbPalette[1]))+
      theme_void()
    return(p)
  }
  
  nodesNoParents <- groupNames[!groupNames%in%edgeList$to]
  vertices <- rbind(vertices,rep(T,dim(vertices)[2]))
  vertices[dim(vertices)[1],c("name","Selected")] <- c("origin",T)
  vertices$visible <- T
  vertices$visible[vertices$name=="origin"] <- F
  edgeList <- rbind(edgeList,data.frame(from=rep("origin",length(nodesNoParents)),to=as.character(nodesNoParents)))
  edgeList <- edgeList %>% mutate(visible=from!="origin")
  mygraph <- graph_from_data_frame( edgeList ,vertices=vertices)
  lay <- create_layout(mygraph, layout = 'dendrogram', circular=F)
  
  #str(get_edges("short")(lay))
  #lay$y[lay$name=="2"]<-lay$y[lay$name=="1"]
  p<- ggraph(lay) + 
    geom_edge_link(aes(alpha=node1.visible),edge_width=ls,
                   arrow=arrow(ends="last",length=unit(3,"mm"),type="closed"),end_cap = circle(6, 'mm')) +
    scale_edge_alpha_discrete(range=c(0,1),guide=F)+
    geom_node_label(aes(label=name2,fill=GrpWeight,col=Selected),
                    show.legend=showFillLegend, fontface = "bold",size=nodeSize)+
    #scale_fill_gradient(low="white",high=cbPalette[4])+
    scale_fill_gradientn(colors=c("white",cols),
                         #breaks=c(0,1,NA),
                         limits=c(0,NA),
                         values = scales::rescale(
                           c(0,min(lay$GrpWeight[lay$GrpWeight>0],10^-6), #white for 0
                             min(10^-6,lay$GrpWeight[lay$GrpWeight>0]),1-1e-6,  #red for <1
                             1-1e-6,1+1e-6,
                             1+1e-6, max(lay$GrpWeight))),
                         name="Group weight")+
    #geom_node_label(aes(label=name2,fill=factor(Selected)), fontface = "bold")+
    scale_color_manual(values=c("FALSE"=cbPalette[1],"TRUE"="black"))+
    coord_cartesian(ylim=c(min(lay$y)-0.5,max(lay$y[lay$name!="origin"])+0.5),xlim=c(min(lay$x)-0.5,max(lay$x)+0.5))+
    #geom_node_text(aes(label=GrpWeight, filter=leaf) , angle=90 , hjust=1, nudge_y = -0.2) +
    theme_void()
  
  return(p)
}

#Visualise grouping weights in CV folds----
visualiseGroupingweights <- function(dfGrps,GroupingNames,hist=F,boxplot=T,jitter=T,ps=1.5,width=0.5){
  #Plot cross-validated grouping weights
  #
  #Input:
  #-dfGrps: dataframe with the following variables:
  #   -Grouping: factor with grouping names
  #   -Grouping.weight: groupings weight of each grouping
  #   -Fold: number indicating which fold in the cross-validation is used
  #   -hist: if hist=T, a histogram is plotted
  #-GroupingNames: vector with grouping names in order levels(dfGrps)
  #
  #Output:
  #-p1: plot in ggplot object
  
  #change grouping names to those in GroupingNames if given
  if(!missing(GroupingNames) && length(GroupingNames)>0){
    if(is.factor(dfGrps$Grouping)){
      dfGrps$Grouping <- factor(dfGrps$Grouping,levels = levels(dfGrps$Grouping)[levels(dfGrps$Grouping)%in%unique(dfGrps$Grouping)],
                                 labels = GroupingNames)
    }else{
      dfGrps$Grouping <- factor(GroupingNames[dfGrps$Grouping],levels=GroupingNames,labels=GroupingNames)
    }
  }else{
    dfGrps$Grouping <- as.factor(dfGrps$Grouping)
  }
  
  if(hist){
    p1 <- ggplot(dfGrps)+
      aes(x=Grouping.weight)+
      geom_histogram(aes(fill=Grouping),position = "identity",bins=20,alpha=0.6)+
      scale_fill_discrete(name="Grouping")+
      labs(x="Grouping weight")+
      theme_bw()+
      theme(axis.text.x=element_text(size=12),
            axis.text.y=element_text(size=12),
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14),
            legend.text=element_text(size=12),
            legend.title=element_text(size=14))
    
  }else if(!(boxplot|jitter)){
    SummdfGrps <- dfGrps %>% group_by(Grouping) %>% summarise("meanWeight"=mean(Grouping.weight),
                                                               "q25Weight"=quantile(Grouping.weight,0.25),
                                                               "q75Weight"=quantile(Grouping.weight,0.75),
                                                               "maxWeight"=max(Grouping.weight),
                                                               "minWeight"=min(Grouping.weight)) %>% ungroup()
    
    p1<-ggplot(SummdfGrps)+
      geom_linerange(aes(x=Grouping,ymax=q75Weight,ymin=q25Weight,col=Grouping),linetype="dashed")+
      geom_point(aes(x=Grouping,y=minWeight,col=Grouping),shape=2)+
      geom_point(aes(x=Grouping,y=meanWeight,col=Grouping),shape=1)+
      geom_point(aes(x=Grouping,y=maxWeight,col=Grouping),shape=6)+
      scale_color_discrete(name="Grouping")+
      labs(y="Grouping weight",x="Grouping")+
      theme_bw()+
      theme(axis.text.x=element_text(size=12),
            axis.text.y=element_text(size=12),
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14),
            legend.text=element_text(size=12),
            legend.title=element_text(size=14))
  }else{
    dfGrps2 <- dfGrps %>% group_by(Fold) %>% distinct(Grouping, .keep_all=T) %>% ungroup()
    p1 <- ggplot(dfGrps2)+
      scale_color_discrete(name="Grouping")+
      labs(y="Grouping weight",x="Grouping")+
      theme_bw()+
      theme(axis.text.x=element_text(size=12),
            axis.text.y=element_text(size=12),
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14),
            legend.text=element_text(size=12),
            legend.title=element_text(size=14))
    if(boxplot){
      p1 <- p1+geom_boxplot(aes(x=Grouping,y=Grouping.weight,col=Grouping),outlier.shape = NA,width=width)
    } 
    if(jitter){
      p1<- p1+geom_jitter(aes(x=Grouping,col=Grouping,y=Grouping.weight),
                          size=ps,alpha=0.6,height=0,width=width/2)
    } 
  }
  
  
  
  return(p1)
}

#Visualise group weights in CV folds----
visualiseGroupweights <- function(dfGrps,Grouping,grouping.grouplvl,values,widthBoxplot=0.05,boxplot=T,jitter=T,ps=1.5,ls=1){
  #Plot cross-validated group weights for one grouping
  #
  #Input:
  #-dfGrps: dataframe with the following variables:
  #   -Group: factor with group names
  #   -Group.weight: group weight of each group
  #   -Fold: number indicating which fold in the cross-validation is used
  #-Grouping: list of G elements containing covariate indices for each group
  #-grouping.grouplvl (optional): groups of groups, e.g. in a hierarchical structure. 
  #-values (optional): values of continuous co-data. If given, group weights are plotted against these value
  #
  #Output:
  #-p1: plot in ggplot object
  
  #change group names to those in Grouping if given
  if(!missing(Grouping) && length(names(Grouping))>0){
    if(is.factor(dfGrps$Group)){
      dfGrps$Group <- factor(dfGrps$Group,levels = levels(dfGrps$Group)[levels(dfGrps$Group)%in%unique(dfGrps$Group)],
                             labels = names(Grouping))
    }else{
      dfGrps$Group <- factor(names(Grouping)[dfGrps$Group])
    }
  }
  
  if(missing(values) || length(values)==0){
    SummdfGrps <- dfGrps %>% group_by(Group) %>% summarise("meanWeight"=mean(Group.weight),
                                                           "q25Weight"=quantile(Group.weight,0.25),
                                                           "q75Weight"=quantile(Group.weight,0.75),
                                                           "maxWeight"=max(Group.weight),
                                                           "minWeight"=min(Group.weight)) %>% ungroup()
    
    if(!(boxplot|jitter)){
      p1<-ggplot(SummdfGrps)+
        geom_hline(yintercept=1,col="grey",linetype="dashed",size=ls,alpha=0.6)+
        geom_linerange(aes(x=Group,ymax=q75Weight,ymin=q25Weight,col=Group),linetype="dashed")+
        geom_point(aes(x=Group,y=minWeight,col=Group),shape=2)+
        geom_point(aes(x=Group,y=meanWeight,col=Group),shape=1)+
        geom_point(aes(x=Group,y=maxWeight,col=Group),shape=6)+
        labs(y="Prior variance weight")+
        theme_bw()+
        theme(axis.text.x=element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.x=element_text(size=14),
              axis.title.y=element_text(size=14),
              legend.text=element_text(size=12),
              legend.title=element_text(size=14))
    }else{
      p1<-ggplot(dfGrps)+
        geom_hline(yintercept=1,col="grey",linetype="dashed",size=ls,alpha=0.6)+
        labs(y="Prior variance weight")+
        theme_bw()+
        theme(axis.text.x=element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.x=element_text(size=14),
              axis.title.y=element_text(size=14),
              legend.text=element_text(size=12),
              legend.title=element_text(size=14))
      if(boxplot){
        p1<- p1+geom_boxplot(aes(x=Group,y=Group.weight,col=Group),outlier.shape = NA,show.legend=F)
      } 
      if(jitter){
        p1<- p1+geom_jitter(aes(x=Group,y=Group.weight,col=Group),show.legend=F,
                            size=ps,alpha=0.6,height=0,width=widthBoxplot/2)
      } 
    }
    return(p1)
  }else{
    getIndBetas<-function(p,ind){
      indBetas <- lapply(1:p,function(x){
        return(unlist(sapply(1:length(ind),function(grp){
          if(any(x==ind[[grp]])) return(grp)
          else return(NULL)
        })))
      }) #for each beta 1,..,p: group numbers
      return(indBetas)
    }
    
    GetWeight <- function(p,indBetas,groupweights,groupnumber){
      averageweight <- unlist(sapply(1:p,function(x){
        grps <- indBetas[[x]]
        size <- length(grps)
        return(mean(groupweights[groupnumber%in%grps]))
      }))
      return(averageweight)
    }
    
    #find leaves (groups that have no children)
    if(all(sapply(obtainHierarchy(Grouping),length)==1)) leaves <- 1:length(Grouping) 
    else leaves <- which(sapply(obtainHierarchy(grouping.grouplvl),length)==1) 
    
    IndBetas<-getIndBetas(p=length(values),Grouping)
    avrg <- sapply(Grouping[leaves],function(x){mean(values[x])}) #average of continuous co-data in leaf groups
    medians <- sapply(Grouping[leaves],function(x){median(values[x])}) #median of continuous co-data in leaf groups
    ord <- order(avrg) 
    leaves<-leaves[ord]
    newLeaf <- as.factor(sapply(IndBetas,function(x){
      which(leaves%in%x)
    }))
    originalLeaf <- as.factor(sapply(IndBetas,function(x){
      leaves[which(leaves%in%x)]
    }))
    originalLeaf <- factor(originalLeaf,levels=leaves,labels=leaves) #order labels
    
    dfBeta <- data.frame()
    for(l in unique(dfGrps$Fold)){
      dfBeta2 <- data.frame("index"=1:length(values),
                            "Weight"=GetWeight(length(values),indBetas=IndBetas,
                                               groupweights=dfGrps$Group.weight[dfGrps$Fold==l],
                                               groupnumber=dfGrps$Group[dfGrps$Fold==l]),
                            "Continuous.values"=values,
                            "newLeafGrp"=newLeaf,
                            "originalLeafGrp"=originalLeaf
      )
      dfBeta2$Fold <- l
      dfBeta2$AverageGroupValue <- avrg[ord[dfBeta2$originalLeafGrp]]
      dfBeta2$MedianGroupValue <- medians[ord[dfBeta2$originalLeafGrp]]
      dfBeta <- rbind(dfBeta,dfBeta2)
    }
    dfBeta$Fold <- as.factor(dfBeta$Fold)
    
    SummdfBeta <- dfBeta %>% group_by(newLeafGrp) %>% summarise("meanValue"=mean(Continuous.values),
                                                                "minValue"=min(Continuous.values),
                                                                "maxValue"=max(Continuous.values),
                                                                "q50Value"=quantile(Continuous.values,0.5),
                                                                "meanWeight"=mean(Weight),
                                                                "q25Weight"=quantile(Weight,0.25),
                                                                "q75Weight"=quantile(Weight,0.75),
                                                                "maxWeight"=max(Weight),
                                                                "minWeight"=min(Weight),
                                                                "q50Weight"=quantile(Weight,0.5),
                                                                "originalLeafGrp"=originalLeafGrp[1]) %>% ungroup()
    
    if(any(is.na(SummdfBeta$meanValue))){ #missing data group
      missingGroup <- which(is.na(SummdfBeta$meanValue))
      SummdfBeta[missingGroup,"minValue"]<-min(SummdfBeta$minValue,na.rm=T)
      SummdfBeta[missingGroup,"maxValue"]<-max(SummdfBeta$maxValue,na.rm=T)
      SummdfBeta[missingGroup,"meanValue"]<-mean(c(SummdfBeta[[missingGroup,"minValue"]],SummdfBeta[[missingGroup,"maxValue"]]))
      levels(SummdfBeta$originalLeafGrp)[missingGroup] <- paste(levels(SummdfBeta$originalLeafGrp)[missingGroup]," (missing data group)",sep="")
    }
    
    if(!(boxplot|jitter)){
      p1<-ggplot(SummdfBeta)+
        geom_hline(yintercept=1,col="grey",linetype="dashed",size=ls,alpha=0.6)+
        geom_linerange(aes(x=q50Value,ymax=q75Weight,ymin=q25Weight,col=originalLeafGrp),linetype="dashed")+
        geom_point(aes(x=q50Value,y=minWeight,col=originalLeafGrp),shape=2)+
        geom_point(aes(x=q50Value,y=meanWeight,col=originalLeafGrp),shape=1)+
        geom_point(aes(x=q50Value,y=maxWeight,col=originalLeafGrp),shape=6)+
        geom_segment(aes(x=minValue,xend=maxValue,y=q50Weight,yend=q50Weight,col=originalLeafGrp))+
        scale_color_discrete(name="Group")+
        labs(x="Continuous co-data value",y="Prior variance weight")+
        theme_bw()+
        theme(axis.text.x=element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.x=element_text(size=14),
              axis.title.y=element_text(size=14),
              legend.text=element_text(size=12),
              legend.title=element_text(size=14))
    }else{
      dfBeta2 <- dfBeta %>% group_by(Fold) %>% distinct(originalLeafGrp, .keep_all=T) %>% ungroup()
      
      p1<-ggplot(dfBeta2)+
        geom_hline(yintercept=1,col="grey",linetype="dashed",size=ls,alpha=0.6)+
        geom_segment(data=SummdfBeta,aes(x=minValue,xend=maxValue,y=q50Weight,yend=q50Weight,col=originalLeafGrp),size=ls)+
        labs(x="Continuous co-data value",y="Prior variance weight")+
        scale_color_discrete(name="Group")+
        theme_bw()+
        theme(axis.text.x=element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.x=element_text(size=14),
              axis.title.y=element_text(size=14),
              legend.text=element_text(size=12),
              legend.title=element_text(size=14))
      if(boxplot){
        p1<- p1 +geom_boxplot(aes(x=MedianGroupValue,col=originalLeafGrp,y=Weight),width=widthBoxplot,
                              varwidth=F,position="identity",outlier.shape = NA)
      } 
      if(jitter){
        p1<- p1+geom_jitter(data=dfBeta2,aes(x=MedianGroupValue,col=originalLeafGrp,y=Weight),
                            size=ps,alpha=0.6,height=0,width=widthBoxplot/4)
      }
      
    }
    return(p1)
  }
}

#get weights of leaf groups in continuous codata
.getWeightLeaves <- function(dfGrps,Grouping,grouping.grouplvl,values){
  #Return data frame with group weights of the leaf groups in a continuous co-data,
  #instead of the group weights per group in the hierarchical tree
  #
  #Input:
  #-dfGrps: dataframe with the following variables:
  #   -Group: factor with group names
  #   -Group.weight: group weight of each group
  #   -Fold: number indicating which fold in the cross-validation is used
  #-Grouping: list of G elements containing covariate indices for each group
  #-grouping.grouplvl (optional): groups of groups, e.g. in a hierarchical structure. 
  #-values (optional): values of continuous co-data. If given, group weights are plotted against these value
  #
  #Output:
  #-df: dataframe similar to dfGrps but with group weights for leaf groups
  
  #change group names to those in Grouping if given
  if(!missing(Grouping) && length(names(Grouping))>0){
    if(is.factor(dfGrps$Group)){
      dfGrps$Group <- factor(dfGrps$Group,levels = levels(dfGrps$Group)[levels(dfGrps$Group)%in%unique(dfGrps$Group)],
                             labels = names(Grouping))
    }else{
      dfGrps$Group <- factor(names(Grouping)[dfGrps$Group])
    }
  }
  
  
    getIndBetas<-function(p,ind){
      indBetas <- lapply(1:p,function(x){
        return(unlist(sapply(1:length(ind),function(grp){
          if(any(x==ind[[grp]])) return(grp)
          else return(NULL)
        })))
      }) #for each beta 1,..,p: group numbers
      return(indBetas)
    }
    
    GetWeight <- function(p,indBetas,groupweights,groupnumber){
      averageweight <- unlist(sapply(1:p,function(x){
        grps <- indBetas[[x]]
        size <- length(grps)
        return(mean(groupweights[groupnumber%in%grps]))
      }))
      return(averageweight)
    }
    
    #find leaves (groups that have no children)
    if(all(sapply(obtainHierarchy(Grouping),length)==1)) leaves <- 1:length(Grouping) 
    else leaves <- which(sapply(obtainHierarchy(grouping.grouplvl),length)==1) 
    
    IndBetas<-getIndBetas(p=length(values),Grouping)
    avrg <- sapply(Grouping[leaves],function(x){mean(values[x])}) #average of continuous co-data in leaf groups
    ord <- order(avrg) 
    leaves<-leaves[ord]
    newLeaf <- as.factor(sapply(IndBetas,function(x){
      which(leaves%in%x)
    }))
    originalLeaf <- as.factor(sapply(IndBetas,function(x){
      leaves[which(leaves%in%x)]
    }))
    originalLeaf <- factor(originalLeaf,levels=leaves,labels=leaves) #order labels
    
    dfBeta <- data.frame()
    for(l in unique(dfGrps$Fold)){
      dfBeta2 <- data.frame("index"=1:length(values),
                            "Weight"=GetWeight(length(values),indBetas=IndBetas,
                                               groupweights=dfGrps$Group.weight[dfGrps$Fold==l],
                                               groupnumber=dfGrps$Group[dfGrps$Fold==l]),
                            "Continuous.values"=values,
                            "newLeafGrp"=newLeaf,
                            "originalLeafGrp"=originalLeaf
      )
      dfBeta2$Fold <- l
      dfBeta <- rbind(dfBeta,dfBeta2)
    }
    dfBeta$Fold <- as.factor(dfBeta$Fold)

    #data frame with the group weight of all betas in the leaf group
    dfGrps2 <- dfBeta %>% group_by(Fold) %>% distinct(originalLeafGrp, .keep_all = TRUE) %>% 
      select(c(originalLeafGrp,Weight,Fold,newLeafGrp)) %>% ungroup()
    dfGrps2$Group <- dfGrps2$originalLeafGrp
    dfGrps2$Grouping<-dfGrps$Grouping[1]
    dfGrps2$Group.weight <- dfGrps2$Weight
    dfGrps2$Method <- dfGrps$Method[1]
    dfGrps2 <- dfGrps2 %>% group_by(Fold) %>% mutate(Grouping.weight=dfGrps$Grouping.weight[which(dfGrps$Fold==unique(Fold))[1]],
                                                     Taugr=dfGrps$Taugr[which(dfGrps$Fold==unique(Fold))[1]],
                                                     Tauglm=dfGrps$Tauglm[which(dfGrps$Fold==unique(Fold))[1]]) %>% ungroup()
    if(!all(names(dfGrps)%in%names(dfGrps2))) warning("Check here, dfGrps has some variables that are not in new data frame")
    dfGrps2 <- dfGrps2 %>% select(names(dfGrps))
    return(dfGrps2)
}