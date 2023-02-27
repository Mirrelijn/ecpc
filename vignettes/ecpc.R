## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo=TRUE
)

## ----setup, eval=FALSE--------------------------------------------------------
#  install.packages("ecpc")

## ----load---------------------------------------------------------------------
library(ecpc)

## ----simdata------------------------------------------------------------------
set.seed(1)
p <- 300 #number of covariates
n <- 100 #sample size training data set
n2 <-100 #sample size test data set
beta <- rnorm(p, mean=0, sd=0.1) #simulate effects
X <- matrix(rnorm(n*p, mean=0, sd=1), n, p) #observed training data
Y <- rnorm(n, mean = X%*%beta, sd=1) #response training data
X2 <- matrix(rnorm(n2*p, mean=0, sd=1), n, p) #observed test data
Y2 <- rnorm(n2, mean = X2%*%beta, sd=1) #response test data

## ----simcodata----------------------------------------------------------------
Z1 <- abs(beta) #informative co-data
Z2 <- rnorm(p, mean=0, sd=1) #random, non-informative co-data
Z <- cbind(Z1, Z2) #(px2)-dimensional co-data matrix 

## -----------------------------------------------------------------------------
fit <- ecpc(Y, X, Z=list(Z), X2=X2, Y2=Y2, postselection=FALSE)

## -----------------------------------------------------------------------------
print(fit$MSEecpc)

## -----------------------------------------------------------------------------
print(fit$MSEridge)

## -----------------------------------------------------------------------------
print(fit)
summary(fit)

## -----------------------------------------------------------------------------
plot(fit, show="coefficients")
plot(fit, show="priorweights", Z=list(Z))

## -----------------------------------------------------------------------------
new_penalties <- penalties(fit, tauglobal = fit$tauglobal * 2, Z=list(Z)) 
new_coefficients <- coef(fit, penalties=new_penalties, X=X, Y=Y) 

## -----------------------------------------------------------------------------
new_penalties2 <- penalties(tauglobal = 1, sigmahat = 1, gamma = c(1,1),
                            w = 1, Z=list(Z)) 
new_coefficients2 <- coef.ecpc(penalties=new_penalties2, X=X, Y=Y) 

## -----------------------------------------------------------------------------
Z1.s <- createZforSplines(values=Z1, G=20, bdeg=3) 
S1.Z1 <- createS(orderPen=2, G=20) 
Z2.s <- createZforSplines(values=Z2, G=30, bdeg=3) 
S1.Z2 <- createS(orderPen=2, G=30)

## -----------------------------------------------------------------------------
Z.all <- list(Z1=Z1.s, Z2=Z2.s)
paraPen.all <- list(Z1=list(S1=S1.Z1), Z2=list(S1=S1.Z2))

## -----------------------------------------------------------------------------
fit.gam <- ecpc(Y, X, Z = Z.all, paraPen = paraPen.all,
                intrcpt.bam=TRUE, X2=X2, Y2=Y2, postselection=FALSE)
fit.gam$MSEecpc


## -----------------------------------------------------------------------------
plot(fit.gam, show="priorweights", Z=Z.all) 
values <- list(Z1, Z2)
plot(fit.gam, show="priorweights", Z=Z.all, values = values) 

## ---- eval=FALSE--------------------------------------------------------------
#  codataNO <- attributes(fit.gam$gamma)$codataSource
#  i <- 2 #1 for informative, 2 for non-informative
#  sk <- as.vector(Z.all[[i]]%*%fit.gam$gamma[codataNO==i])*fit.gam$tauglobal
#  par(mfrow=c(1,1))
#  plot(Z[,i],sk)

## -----------------------------------------------------------------------------
Con.Z1 <- createCon(G=20, shape="positive+monotone.i") 
Con.Z2 <- createCon(G=30, shape="convex") 
paraCon <- list(Z1=Con.Z1, Z2=Con.Z2)

## -----------------------------------------------------------------------------
fit.scam <- ecpc(Y, X, Z = Z.all, paraPen = paraPen.all, 
                 paraCon = paraCon, X2=X2, Y2=Y2, postselection=FALSE)
fit.scam$MSEecpc
plot(fit.scam, show="priorweights", Z=Z.all, values=values)

## -----------------------------------------------------------------------------
Ypred <- predict(fit, X2=X2)

## ---- eval=FALSE--------------------------------------------------------------
#  if(!requireNamespace("squeezy")) install.packages("squeezy")

## -----------------------------------------------------------------------------
library("squeezy")
fit.EN <- squeezy(Y, X, alpha=0.3, X2=X2, Y2=Y2, lambdas=fit$penalties)
summary(fit.EN$lambdapApprox) #transformed elastic net penalties
summary(fit.EN$betaApprox) #fitted elastic net regression coefficients

## -----------------------------------------------------------------------------
maxsel= c(5,10)
sparseModels <- postSelect(fit, X=X, Y=Y, maxsel=maxsel)

