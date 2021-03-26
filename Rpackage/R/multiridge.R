#The functions below are contained in the multiridge package, available on github
#Author: Mark van de Wiel

#FITTING: efficient IWLS algorithm for logistic and linear (only one iteration)
.IWLSridge <- function(XXT,Y,X1=NULL,intercept=TRUE,
                       frac1=NULL,eps=1e-7,maxItr=25,trace=FALSE, 
                       model=NULL, E0=NULL){
  #Input:
  #XXT:  sample cross-product (nxn matrix) from penalized variables [computed by SigmaFromBlocks]
  #Y: nx1 response
  #X1: optional nxp1 matrix (p1 < n) with unpenalized covariates
  #intercept: use intercept? If yes, intercept is either estimated or set equal to logit(frac1)
  #frac1: fraction of 1's in the population; if non-null the intercept is estimated as logit(frac1); otherwise intercept is estimated from the data along with all other parameters
  #eps: numerical bound for convergence
  #maxItr: maximum number of iterations used in IWLS
  #trace: should output be traced?
  
  #Output: List object containing:Ypred=Ypred,convergence=convergence,nIt=it,Hres=Hres, linearized=linearized, unpen=unpen,eta0=eta0))
  #etas: linear predictors, vector (n)
  #Ypred: predicted responses, vector (n)
  #convergence (T/F): T if IWLS has converged under threshold eps, F otherwise
  #nIt: number of iterations needed for convergence
  #Hres: list containing matrices required for predictions on new samples
  #linearized: linearized response vector (n) of last IWLS iteration, required for predictions
  #unpen: does the model include unpenalized parameters, including the intercept?
  #eta0: scalar. Offset = logit(frac1), if frac1 is specified. The same for all samples. Required for predictions.
  
  
  #intercept<-FALSE
  #frac1 <- NULL
  #X1 <- matrix(rep(1,n),ncol=1)
  #XXblocks <- list(XXTmir);penalties <- c(mirlambda)  #CIN3=1, normal=1
  #XXblocks <- list(XXTmir,XXTmeth);penalties <- c(mirlambda,methlambda)  #CIN3=1, normal=1
  #maxItr <- 10
  #eps=1e-5
  
  #in terms of hat matrix
  #note: eta is lin pred without intercept; etatot incl intercept
  if(is.null(model)) model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
  
  n <- length(Y)
  if(intercept & is.null(frac1)) X1 <- cbind(matrix(rep(1,n),ncol=1),X1) #if frac1 is not given intercept is estimated alongside all other parameters
  if(intercept & !is.null(frac1)) eta0 <- log((frac1)/(1-frac1)) else eta0 <-  0 #if frac1 is given it is used to estimate the intercept
  
  if(is.null(E0)) eta <- rep(0,n) else eta <- E0
  
  etatot <- eta0 + eta
  it <- 1; nrm <- Inf
  
  if(model=="logistic") Ypred<-1/(1+exp(-etatot))
  if(model=="linear") {Ypred <- etatot; maxItr<-1}
  while(nrm>eps & it<=maxItr){
    etaold <- eta
    Ypredold <- Ypred
    if(model=="logistic") WV <-Ypred*(1-Ypred) #Var(Y)
    if(model=="linear") WV <- rep(1,n)
    #Winv <- diag(1/(Ypred*(1-Ypred)))
    #XtXDinv <- solve(t(X) %*% W %*% X + 2*diag(lambda,p)) #inverting pxp matrix
    znew <-  matrix(Y - Ypred,ncol=1)
    # inv <- solve(Winv + XXT)
    # Hmat <- XXT -  XXT %*% inv %*% XXT
    if(is.null(X1)) {Hres <- .Hpen(WV,XXT);Hmat <- Hres$Hmat; unpen <- FALSE} else {
      Hres <- .Hunpen(WV,XXT,X1); Hmat <- Hres$Hmat; unpen <- TRUE
    }
    # if(is.null(X1)) Hmat <- .Hpen(WV,XXT) else {
    #    Hmat <- .Hunpen(WV,XXT,X1);Hres <- NULL
    #  }
    #eta <- Hmat %*% znew + XXT %*% inv %*% etaold
    eta <- Hmat %*% (znew + diag(WV) %*% etaold)
    etatot <- eta0 + eta
    if(model=="logistic") Ypred<-.minmaxf(1/(1+exp(-etatot))) #predicted probability
    #if(model=="logistic") Ypred<-1/(1+exp(-etatot)) #predicted probability
    if(model=="linear") Ypred <- etatot
    #if(trace) {print(Ypred); print(WV)}
    if(trace) {print(Ypred)}
    it <- it+1
    nrm <- max(abs(Ypred-Ypredold))
    if(trace) print(nrm)
    #alt verified for mir data
    # eta2 <- datamir %*% solve(t(datamir) %*% W %*% datamir + diag(rep(penalties[1],ncol(datamir)))) %*% t(datamir) %*% z
  }
  linearized <- znew + diag(WV) %*% etaold
  if(trace) print(it)
  convergence <- nrm<eps
  return(list(etas=eta,Ypred=Ypred,convergence=convergence,nIt=it-1,
              Hres=Hres, linearized=linearized, unpen=unpen,
              intercept=intercept, eta0=eta0)) #KW, MW components of Hres needed for pred, and so is the linearized response
}

#IWLS for Cox ridge
.IWLSCoxridge <- function(XXT,Y,X1=NULL,intercept=FALSE,eps=1e-7,maxItr=25,trace=FALSE, E0=NULL){
  #Input:
  #XXT:  sample cross-product (nxn matrix) from penalized variables [computed by SigmaFromBlocks]
  #Y: survival response
  #X1: optional nxp1 matrix (p1 < n) with unpenalized covariates
  #intercept: use intercept? If yes, intercept is either estimated or set equal to logit(frac1)
  #frac1: fraction of 1's in the population; if non-null the intercept is estimated as logit(frac1); otherwise intercept is estimated from the data along with all other parameters
  #eps: numerical bound for convergence
  #maxItr: maximum number of iterations used in IWLS
  #trace: should output be traced?
  
  #Output: List object containing:Ypred=Ypred,convergence=convergence,nIt=it,Hres=Hres, linearized=linearized, unpen=unpen,eta0=eta0))
  #etas: linear predictors, vector (n)
  #Ypred: predicted responses, vector (n)
  #convergence (T/F): T if IWLS has converged under threshold eps, F otherwise
  #nIt: number of iterations needed for convergence
  #Hres: list containing matrices required for predictions on new samples
  #linearized: linearized response vector (n) of last IWLS iteration, required for predictions
  #unpen: does the model include unpenalized parameters, including the intercept?
  #eta0: scalar. Offset = logit(frac1), if frac1 is specified. The same for all samples. Required for predictions.
  
  
  #intercept<-FALSE
  #frac1 <- NULL
  #X1 <- matrix(rep(1,n),ncol=1)
  #XXblocks <- list(XXTmir);penalties <- c(mirlambda)  #CIN3=1, normal=1
  #XXblocks <- list(XXTmir,XXTmeth);penalties <- c(mirlambda,methlambda)  #CIN3=1, normal=1
  #maxItr <- 10
  #eps=1e-5
  
  #in terms of hat matrix
  #note: eta is lin pred without intercept; etatot incl intercept
  #XXT = SigmaFromBlocks(XXomics[2],penalties=lambdas[2]);Y=resp;X1=NULL;intercept=FALSE;frac1=NULL;eps=1e-7;maxItr=100;trace=TRUE
  
  n <- dim(Y)[1]
  eta0 <- 0
  if(is.null(E0)) eta <- rep(0,n) else eta <- E0
  
  it <- 1; nrm <- Inf
  etatot <- eta
  H0 <- .breslow(Y,etatot)[,2]
  Ypred<- as.numeric(H0 * exp(etatot))
  Yev <- Y[,2]
  while(nrm>eps & it<=maxItr){
    etaold <- eta
    Ypredold <- Ypred
    WV <- Ypred + 10^(-10)
    znew <-  matrix(Yev - Ypred,ncol=1)
    if(is.null(X1)) {Hres <- .Hpen(WV,XXT);Hmat <- Hres$Hmat; unpen <- FALSE} else {
      Hres <- .Hunpen(WV,XXT,X1); Hmat <- Hres$Hmat; unpen <- TRUE
    }
    eta <- Hmat %*% (znew + diag(WV) %*% etaold)
    etatot <- eta
    H0 <- .breslow(Y,etatot)[,2]
    Ypred<- as.numeric(H0 * exp(etatot))
    if(trace) {print(etatot)}
    #if(trace) {print(WV)}
    it <- it+1
    nrm <- max(abs(Ypred-Ypredold))
    if(trace) print(nrm)
  }
  linearized <- znew + diag(WV) %*% etaold
  if(trace) print(it)
  convergence <- nrm<eps
  return(list(etas=eta,Ypred=Ypred,convergence=convergence,nIt=it-1,
              Hres=Hres, linearized=linearized, unpen=unpen,
              intercept=intercept,eta0=eta0)) #KW, MW components of Hres needed for pred, and so is the linearized response
}

#Computation of weighted Hat matrix, penalized variables only; needed in IWLS functions
.Hpen <- function(WV,XXT){ #WV:weigths as vector (n); XXT: (penalized) sample cross-product (nxn)
  n <-length(WV)
  Winv <- diag(1/WV)
  inv <- solve(Winv + XXT)
  
  Mmat <- diag(n) - inv %*% XXT
  Hmat <- XXT %*% Mmat
  return(list(Hmat=Hmat,Mmat=Mmat))
}

#Computation of weighted Hat matrix, accounting for unpenalized variables; needed in IWLS functions
.Hunpen <- function(WV,XXT,X1){
  #WV:weigths as vector (n)
  #XXT: sample cross-product (nxn) penalized variables
  #X1 (p1xn) design matrix unpenalized variables
  n <- length(WV)
  WsqrtV <- sqrt(WV)
  WinvsqrtV <- 1/WsqrtV
  X1W <- WsqrtV * X1
  X1aux <- solve(t(X1W) %*% X1W) %*% t(X1W)
  X1Wproj <- X1W %*% X1aux
  P1W <- diag(n) - X1Wproj
  GammaW <- t(t(WsqrtV * XXT) * WsqrtV) #faster
  P1GammaW <- P1W %*% GammaW

  invforH2<-try(solve(diag(n) + P1GammaW))
  if(class(invforH2)[1]=="try-error"){
    svdXXT <- svd(diag(n) + P1GammaW)
    svdd <- svdXXT$d
    #reci <- 1/svdd[1:n]
    reci <- c(1/svdd[1:n-1],0)
    invforH2 <- svdXXT$v %*% (reci * t(svdXXT$u))
  }
  
  #H2 <- Winvsqrt %*% GammaW %*% (diag(n) -invforH2 %*% P1GammaW) %*% P1W %*% Winvsqrt
  MW <- WsqrtV * t(t((diag(n) - invforH2 %*% P1GammaW) %*% P1W) * WinvsqrtV)
  H20 <- GammaW %*% (WinvsqrtV * MW)
  H2 <- t(t(WinvsqrtV * H20)) #faster
  #Hboth <- Winvsqrt %*% X1Wproj %*% (diag(n) - Wsqrt %*% H2 %*% Wsqrt) %*% Winvsqrt + H2
  KW <- t(t(X1aux %*% (diag(n) - t(t(WsqrtV * H2) * WsqrtV))) * WinvsqrtV)
  Hmat <-(WinvsqrtV * X1W) %*% KW  + H2
  return(list(Hmat=Hmat,MW=MW,KW=KW))
}

#auxiliary function; thresholds predictions to avoid numerical problems; in IWLSRidge
.minmaxf <- function(vec,thr=10^(-3)){
  sapply(vec,function(x) min(max(x,thr),1-thr))
}

#NEEDED FOR IWLSCoxRidge
.breslow <- function(Y,lp){ #checked to coincide with basehaz from penalized (which coincides with survfit)
  #Y survival response; lp: linear predictor (length n)
  #Returns hazard and cumulative hazard at all timepoints
  ord <- order(Y[,1])
  invord <- order(ord) #to put back in original order
  Ys <- Y[ord,]
  di <- Ys[,2] #event
  ti <- Ys[,1] #time
  lp <- lp[ord]
  htsort <- di/rev(cumsum(exp(rev(lp))))
  ht <- htsort[invord]
  Ht <- cumsum(htsort)[invord]
  return(data.frame(ht=ht,Ht=Ht))
}

#Predictions for new samples
.predictIWLS <- function(IWLSfit,X1new=NULL, Sigmanew){
  #IWLSfit: list object; output from IWLSridge
  #X1new: design matrix with unpenalized variables, p1 x nnew; should not contain the intercept
  #Sigmanew: (penalized) sample cross-product between new and training samples (nnew x n); may be computed using SigmaFromBlocks function
  #Output: vector (nnew) of linear predictors
  
  unpen <- IWLSfit$unpen
  eta0 <- IWLSfit$eta0
  intercept <- IWLSfit$intercept
  if(unpen){
    nout <- nrow(Sigmanew)
    if(intercept) X1new <- cbind(matrix(rep(1,nout),ncol=1),X1new)
    linearized <- IWLSfit$linearized
    KW <- IWLSfit$Hres$KW
    MW <- IWLSfit$Hres$MW
    if(ncol(X1new) != nrow(KW)){
      print("Error: Unpenalized design matrix should contain same number of columns as Unpenalized design matrix used for training")
      return(NULL)
    } else Hnew <- X1new %*% KW + Sigmanew %*% MW
  } else {
    linearized <- IWLSfit$linearized
    Mmat <- IWLSfit$Hres$Mmat
    Hnew <-Sigmanew %*% Mmat
  }
  pred <- Hnew %*% linearized + eta0
  return(pred)
}