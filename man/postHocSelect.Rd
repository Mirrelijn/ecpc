\name{postHocSelect}
\alias{postHocSelect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Perform posterior selection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
postHocSelect(X, Y, beta, intrcpt = 0, penfctr, postselection = "elnet+dense", maxsel = 30, penalties, model, tauglobal, sigmahat = 1, muhatp = 0, X2 = NaN, Y2 = NaN)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{intrcpt}{
%%     ~~Describe \code{intrcpt} here~~
}
  \item{penfctr}{
%%     ~~Describe \code{penfctr} here~~
}
  \item{postselection}{
%%     ~~Describe \code{postselection} here~~
}
  \item{maxsel}{
%%     ~~Describe \code{maxsel} here~~
}
  \item{penalties}{
%%     ~~Describe \code{penalties} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{tauglobal}{
%%     ~~Describe \code{tauglobal} here~~
}
  \item{sigmahat}{
%%     ~~Describe \code{sigmahat} here~~
}
  \item{muhatp}{
%%     ~~Describe \code{muhatp} here~~
}
  \item{X2}{
%%     ~~Describe \code{X2} here~~
}
  \item{Y2}{
%%     ~~Describe \code{Y2} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Y, beta, intrcpt = 0, penfctr, postselection = "elnet+dense", 
    maxsel = 30, penalties, model, tauglobal, sigmahat = 1, muhatp = 0, 
    X2 = NaN, Y2 = NaN) 
{
    tautrgt <- tauglobal
    n <- dim(X)[1]
    p <- dim(X)[2]
    if (missing(penfctr)) 
        penfctr <- rep(1, p)
    maxsel2 <- pmin(maxsel, p)
    if (any(maxsel2 < 2)) {
        warning("Number of variables to be selected should be at least 2 (out of convenience)")
        print("Maxsel values smaller than 2 are set to 2")
        maxsel2[maxsel2 < 2] <- 2
    }
    nonzeros <- beta != 0
    lambdap <- penalties
    pen <- which(penfctr != 0)
    switch(model, linear = {
        fml <- "gaussian"
        sd_y <- sqrt(var(Y) * (n - 1)/n)[1]
    }, logistic = {
        fml <- "binomial"
        sd_y <- 1
    }, cox = {
        fml <- "cox"
        sd_y <- 1
    })
    if (grepl("elnet", postselection)) {
        output <- lapply(maxsel2, function(x) {
            if (sum(beta[pen] != 0) <= x) {
                betaPost <- beta
                output <- list()
                output$betaPost <- betaPost
                output$whichPost <- which(nonzeros & (1:p) \%in\% 
                  pen)
                output$a0 <- intrcpt
                if (!all(is.nan(X2))) {
                  if (model == "linear") {
                    YpredPost <- X2 \%*\% c(betaPost) + intrcpt
                    MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                  }
                  if (model == "logistic") {
                    YpredPost <- 1/(1 + exp(-(X2 \%*\% c(betaPost) + 
                      intrcpt)))
                    MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                    if (any(is.nan(YpredPost))) {
                      browser()
                    }
                  }
                  else if (model == "cox") {
                    expXbPost <- exp(X \%*\% c(betaPost))
                    h0 <- sapply(1:length(Y[, 1]), function(i) {
                      Y[i, 2]/sum(expXbPost[Y[, 1] >= Y[i, 1]])
                    })
                    H0 <- sapply(Y2[, 1], function(Ti) {
                      sum(h0[Y[, 1] <= Ti])
                    })
                    YpredPost <- H0 * exp(X2 \%*\% c(betaPost))
                    MSEPost <- sum((YpredPost - Y2[, 2])^2)/length(Y2[, 
                      2])
                  }
                  output$MSEPost <- MSEPost
                  output$YpredPost <- YpredPost
                }
                return(output)
            }
            else {
                if (length(intrcpt) == 0 || intrcpt == 0) {
                  intrcpt <- F
                }
                else {
                  intrcpt <- T
                }
                if (length(muhatp) == 1) 
                  muhatp <- rep(muhatp, length(pen))
                if (all(muhatp == 0)) 
                  offset <- rep(0, n)
                else offset <- X[, pen] \%*\% muhatp
                lambdaoverall <- sigmahat/tautrgt
                lam2 <- sigmahat/tautrgt/n * sd_y
                Xacc <- X
                Xacc[, pen] <- as.matrix(X[, pen] \%*\% sparseMatrix(i = 1:length(lambdap), 
                  j = 1:length(lambdap), x = c(1/sqrt(lambdap/lambdaoverall))))
                fsel <- function(alpha, maxselec = x) {
                  if (alpha == 0) 
                    return(p - maxselec)
                  else {
                    if (model == "cox") {
                      glmPost <- glmnet(Xacc[, nonzeros], Y, 
                        alpha = alpha, lambda = lam2/(1 - alpha), 
                        family = fml, offset = offset, standardize = F, 
                        penalty.factor = penfctr[nonzeros], thresh = 10^-10)
                    }
                    else if (model \%in\% c("logistic", "linear")) {
                      glmPost <- glmnet(Xacc[, nonzeros], Y, 
                        alpha = alpha, lambda = lam2/(1 - alpha), 
                        family = fml, offset = offset, standardize = F, 
                        intercept = intrcpt, penalty.factor = penfctr[nonzeros], 
                        thresh = 10^-10)
                    }
                    betaPost <- rep(0, p)
                    betaPost[nonzeros] <- as.vector(glmPost$beta)
                    betaPost[pen] <- c(1/sqrt(lambdap/lambdaoverall)) * 
                      betaPost[pen] + muhatp
                    return(sum(betaPost[pen] != 0) - maxselec)
                  }
                }
                rangeAlpha <- c(0, 10/11)
                ItrAlp <- 1
                while (sign(fsel(0)) == sign(fsel(rangeAlpha[2])) & 
                  ItrAlp <= 50) {
                  rangeAlpha[2] <- rangeAlpha[2] + 0.5 * (1 - 
                    rangeAlpha[2])
                  ItrAlp <- ItrAlp + 1
                }
                alpha <- uniroot(fsel, interval = rangeAlpha, 
                  maxiter = 200, tol = 10^(-10))$root
                if (model == "cox") {
                  glmPost0 <- glmnet(Xacc[, nonzeros], Y, alpha = alpha, 
                    lambda = lam2/(1 - alpha), family = fml, 
                    offset = offset, standardize = F, penalty.factor = penfctr[nonzeros], 
                    thresh = 10^-10)
                }
                else {
                  glmPost0 <- glmnet(Xacc[, nonzeros], Y, alpha = alpha, 
                    lambda = lam2/(1 - alpha), family = fml, 
                    offset = offset, intercept = intrcpt, standardize = F, 
                    penalty.factor = penfctr[nonzeros], thresh = 10^-10)
                }
                betaPost0 <- rep(0, p)
                betaPost0[nonzeros] <- as.vector(glmPost0$beta)
                betaPost0[pen] <- c(1/sqrt(lambdap/sigmahat * 
                  tautrgt)) * betaPost0[pen] + muhatp
                whichPostboth <- betaPost0 != 0
                if (sum(whichPostboth) <= 1) {
                  warning("At least two variables should be selected for glmnet")
                  return(list(betaPost = NULL, whichPost = NULL, 
                    a0 = NULL))
                }
                if (grepl("dense", postselection)) {
                  if (!all(muhatp == 0)) 
                    offset <- X[, whichPostboth & (1:p) \%in\% 
                      pen, drop = F] \%*\% muhatp[whichPostboth[pen], 
                      drop = F]
                  if (grepl("dense2", postselection)) {
                    if (model == "cox") {
                      lambdaGLM <- cv.glmnet(Xacc[, whichPostboth, 
                        drop = F], Y, alpha = 0, family = fml, 
                        offset = offset, standardize = F, penalty.factor = penfctr[whichPostboth], 
                        thresh = 10^-10)
                      lam2 <- lambdaGLM$lambda.min
                    }
                    else {
                      lambdaGLM <- cv.glmnet(Xacc[, whichPostboth, 
                        drop = F], Y, alpha = 0, family = fml, 
                        offset = offset, intercept = intrcpt, 
                        standardize = F, penalty.factor = penfctr[whichPostboth], 
                        thresh = 10^-10)
                      lam2 <- lambdaGLM$lambda.min
                    }
                  }
                  if (model == "cox") {
                    glmPost <- glmnet(Xacc[, whichPostboth, drop = F], 
                      Y, alpha = 0, lambda = lam2, family = fml, 
                      offset = offset, standardize = F, penalty.factor = penfctr[whichPostboth], 
                      thresh = 10^-10)
                  }
                  else {
                    glmPost <- glmnet(Xacc[, whichPostboth, drop = F], 
                      Y, alpha = 0, lambda = lam2, family = fml, 
                      offset = offset, intercept = intrcpt, standardize = F, 
                      penalty.factor = penfctr[whichPostboth], 
                      thresh = 10^-10)
                  }
                  betaPost <- rep(0, p)
                  betaPost[whichPostboth] <- as.vector(glmPost$beta)
                  indPostpen <- whichPostboth[pen]
                  indPostp <- whichPostboth & ((1:p) \%in\% pen)
                  betaPost[indPostp] <- c(1/sqrt(lambdap[indPostpen]/sigmahat * 
                    tautrgt)) * betaPost[indPostp] + muhatp[indPostpen]
                  whichPost <- which(indPostp)
                  output <- list()
                  output$betaPost <- betaPost
                  output$whichPost <- whichPost
                  output$a0 <- glmPost$a0
                }
                else {
                  if (model == "cox") {
                    lambdaGLM <- cv.glmnet(X[, whichPostboth, 
                      drop = F], Y, alpha = 0, family = fml, 
                      standardize = F, penalty.factor = penfctr[whichPostboth])
                    lam2 <- lambdaGLM$lambda.min
                  }
                  else {
                    lambdaGLM <- cv.glmnet(X[, whichPostboth, 
                      drop = F], Y, alpha = 0, family = fml, 
                      standardize = F, penalty.factor = penfctr[whichPostboth])
                    lam2 <- lambdaGLM$lambda.min
                  }
                  if (model == "cox") {
                    glmPost <- glmnet(X[, whichPostboth, drop = F], 
                      Y, alpha = 0, lambda = lam2, family = fml, 
                      offset = offset, standardize = F, penalty.factor = penfctr[whichPostboth], 
                      thresh = 10^-10)
                  }
                  else {
                    glmPost <- glmnet(X[, whichPostboth, drop = F], 
                      Y, alpha = 0, lambda = lam2, family = fml, 
                      intercept = intrcpt, standardize = F, penalty.factor = penfctr[whichPostboth], 
                      thresh = 10^-10)
                  }
                  betaPost <- rep(0, p)
                  betaPost[whichPostboth] <- as.vector(glmPost$beta)
                  indPostpen <- whichPostboth[pen]
                  indPostp <- whichPostboth & ((1:p) \%in\% pen)
                  whichPost <- which(indPostp)
                  output <- list()
                  output$betaPost <- betaPost
                  output$whichPost <- whichPost
                  output$a0 <- glmPost$a0
                }
                if (!all(is.nan(X2))) {
                  if (model == "linear") {
                    YpredPost <- X2 \%*\% c(betaPost) + output$a0
                    MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                  }
                  if (model == "logistic") {
                    YpredPost <- 1/(1 + exp(-(X2 \%*\% c(betaPost) + 
                      output$a0)))
                    MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                    if (any(is.nan(YpredPost))) {
                      browser()
                    }
                  }
                  else if (model == "cox") {
                    expXbPost <- exp(X \%*\% c(betaPost))
                    h0 <- sapply(1:length(Y[, 1]), function(i) {
                      Y[i, 2]/sum(expXbPost[Y[, 1] >= Y[i, 1]])
                    })
                    H0 <- sapply(Y2[, 1], function(Ti) {
                      sum(h0[Y[, 1] <= Ti])
                    })
                    YpredPost <- H0 * exp(X2 \%*\% c(betaPost))
                    MSEPost <- sum((YpredPost - Y2[, 2])^2)/length(Y2[, 
                      2])
                  }
                  output$MSEPost <- MSEPost
                  output$YpredPost <- YpredPost
                }
                return(output)
            }
        })
    }
    else if (grepl("DSS", postselection)) {
        Xacc <- X
        Xacc[, pen] <- as.matrix(Xacc[, pen] \%*\% sparseMatrix(i = 1:length(beta[pen]), 
            j = 1:length(pen), x = c(sqrt(abs(beta[pen])))))
        Ygamma <- Xacc \%*\% beta
        if (grepl("fast", postselection)) {
            glmPost <- glmnet(x = Xacc[, nonzeros], y = Ygamma, 
                alpha = 1, nlambda = 100, family = "gaussian", 
                intercept = intrcpt, standardize = F, penalty.factor = penfctr[nonzeros], 
                thresh = 10^-10)
        }
        else {
            if (length(intrcpt == 0) || intrcpt == 0) {
                intrcpt <- F
            }
            else {
                intrcpt <- T
            }
            fsel <- function(lam1, maxselec = maxsel2) {
                glmPost <- glmnet(x = Xacc[, nonzeros], y = Ygamma, 
                  alpha = 1, lambda = lam1, family = "gaussian", 
                  intercept = intrcpt, standardize = F, penalty.factor = penfctr[nonzeros], 
                  thresh = 10^-10)
                return(sum(as.vector(glmPost$beta) != 0) - maxselec)
            }
        }
        output <- lapply(maxsel2, function(x) {
            if (length(beta[pen] != 0) <= x) {
                betaPost <- beta
                output <- list()
                output$betaPost <- betaPost
                output$whichPost <- which(betaPost != 0)
                output$a0 <- intrcpt
                if (!all(is.nan(X2))) {
                  if (model == "linear") {
                    YpredPost <- X2 \%*\% c(betaPost) + intrcpt
                    MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                  }
                  if (model == "logistic") {
                    YpredPost <- 1/(1 + exp(-(X2 \%*\% c(betaPost) + 
                      intrcpt)))
                    MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                    if (any(is.nan(YpredPost))) {
                      browser()
                    }
                  }
                  else if (model == "cox") {
                    expXbPost <- exp(X \%*\% c(betaPost))
                    h0 <- sapply(1:length(Y[, 1]), function(i) {
                      Y[i, 2]/sum(expXbPost[Y[, 1] >= Y[i, 1]])
                    })
                    H0 <- sapply(Y2[, 1], function(Ti) {
                      sum(h0[Y[, 1] <= Ti])
                    })
                    YpredPost <- H0 * exp(X2 \%*\% c(betaPost))
                    MSEPost <- sum((YpredPost - Y2[, 2])^2)/length(Y2[, 
                      2])
                  }
                  output$MSEPost <- MSEPost
                  output$YpredPost <- YpredPost
                }
                return(output)
            }
            else {
                if (grepl("fast", postselection)) {
                  whlam <- which(glmPost$df <= x)
                  takelam <- rev(whlam)[1]
                  lam <- glmPost$lambda[takelam]
                  betaPost <- rep(0, p)
                  betaPost[nonzeros] <- as.vector(glmPost$beta[, 
                    takelam])
                  betaPost[pen] <- c(sqrt(abs(beta[pen]))) * 
                    betaPost[pen]
                  whichPost <- which(betaPost != 0 & ((1:p) \%in\% 
                    pen))
                  a0 <- glmPost$a0[takelam]
                }
                else {
                  lam1 <- uniroot(fsel, maxselec = x, interval = c(0, 
                    max(lambdap[lambdap < Inf]) * sigmahat/n), 
                    maxiter = 200, tol = 10^(-10))$root
                  glmPost <- glmnet(x = Xacc[, nonzeros], y = Ygamma, 
                    alpha = 1, lambda = lam1, family = "gaussian", 
                    intercept = intrcpt, standardize = F, penalty.factor = penfctr[nonzeros], 
                    thresh = 10^-10)
                  betaPost <- rep(0, p)
                  betaPost[nonzeros] <- as.vector(glmPost$beta)
                  betaPost[pen] <- c(sqrt(abs(beta[pen]))) * 
                    betaPost[pen]
                  whichPost <- which(betaPost != 0 & ((1:p) \%in\% 
                    pen))
                  a0 <- glmPost$a0
                }
                output <- list()
                output$betaPost <- betaPost
                output$whichPost <- whichPost
                output$a0 <- a0
                if (!all(is.nan(X2))) {
                  if (model == "linear") {
                    YpredPost <- X2 \%*\% c(betaPost) + output$a0
                    MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                  }
                  if (model == "logistic") {
                    YpredPost <- 1/(1 + exp(-(X2 \%*\% c(betaPost) + 
                      output$a0)))
                    MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                    if (any(is.nan(YpredPost))) {
                      browser()
                    }
                  }
                  else if (model == "cox") {
                    expXbPost <- exp(X \%*\% c(betaPost))
                    h0 <- sapply(1:length(Y[, 1]), function(i) {
                      Y[i, 2]/sum(expXbPost[Y[, 1] >= Y[i, 1]])
                    })
                    H0 <- sapply(Y2[, 1], function(Ti) {
                      sum(h0[Y[, 1] <= Ti])
                    })
                    YpredPost <- H0 * exp(X2 \%*\% c(betaPost))
                    MSEPost <- sum((YpredPost - Y2[, 2])^2)/length(Y2[, 
                      2])
                  }
                  output$MSEPost <- MSEPost
                  output$YpredPost <- YpredPost
                }
                return(output)
            }
        })
    }
    else if (grepl("BR", postselection)) {
        if (grepl("joint", postselection)) {
            print("Bondell&Reich joint credible set used to select covariates")
            ind <- which(penfctr != 0 & penalties != Inf)
            if (model == "linear") {
                W <- diag(rep(1, n))
            }
            else if (model == "logistic") {
                expminXb <- exp(-X \%*\% c(beta) - intrcpt)
                Pi <- 1/(1 + expminXb)
                W <- diag(c(sqrt(Pi * (1 - Pi))))
            }
            else if (model == "cox") {
                expXb <- exp(X \%*\% c(beta))
                h0 <- sapply(1:length(Y[, 1]), function(i) {
                  Y[i, 2]/sum(expXb[Y[, 1] >= Y[i, 1]])
                })
                H0 <- sapply(Y[, 1], function(Ti) {
                  sum(h0[Y[, 1] <= Ti])
                })
                W <- diag(c(sqrt(H0 * expXb)))
            }
            Sigminhalf <- sqrtm(t(X[, ind]) \%*\% W \%*\% X[, ind] + 
                diag(lambdap[ind]))/sqrt(sigmahat)
            Xstar <- t(t(Sigminhalf) * beta[ind]^2)
            Ystar <- Sigminhalf \%*\% beta[ind]
            glmPost <- glmnet(Xstar, Ystar, alpha = 1, family = "gaussian", 
                standardize = F, intercept = F, thresh = 10^-10)
            output <- lapply(maxsel2, function(x) {
                fsel <- function(lambda, maxselec = x) {
                  fit <- coef.glmnet(glmPost, s = lambda, exact = T, 
                    x = Xstar, y = Ystar)
                  return(sum(fit != 0) - maxselec)
                }
                lambda <- uniroot(fsel, interval = c(0, max(glmPost$lambda)), 
                  maxiter = 200, tol = 10^(-10))$root
                fit <- coef.glmnet(glmPost, s = lambda, exact = T, 
                  x = Xstar, y = Ystar)
                if (grepl("same", postselection)) {
                  if (x == maxsel2[1]) 
                    print("Selected covariates are not refit, unpenalised covariates are kept the same")
                  betaPost <- beta
                  betaPost[ind] <- fit * beta[ind]^2
                  whichPost <- which(betaPost != 0 & penfctr != 
                    0)
                  output <- list()
                  output$betaPost <- betaPost
                  output$whichPost <- whichPost
                  output$a0 <- intrcpt
                  if (!all(is.nan(X2))) {
                    if (model == "linear") {
                      YpredPost <- X2 \%*\% c(betaPost) + intrcpt
                      MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                    }
                    if (model == "logistic") {
                      YpredPost <- 1/(1 + exp(-(X2 \%*\% c(betaPost) + 
                        intrcpt)))
                      MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                      if (any(is.nan(YpredPost))) {
                        browser()
                      }
                    }
                    else if (model == "cox") {
                      expXbPost <- exp(X \%*\% c(betaPost))
                      h0 <- sapply(1:length(Y[, 1]), function(i) {
                        Y[i, 2]/sum(expXbPost[Y[, 1] >= Y[i, 
                          1]])
                      })
                      H0 <- sapply(Y2[, 1], function(Ti) {
                        sum(h0[Y[, 1] <= Ti])
                      })
                      YpredPost <- H0 * exp(X2 \%*\% c(betaPost))
                      MSEPost <- sum((YpredPost - Y2[, 2])^2)/length(Y2[, 
                        2])
                    }
                    output$MSEPost <- MSEPost
                    output$YpredPost <- YpredPost
                  }
                  return(output)
                }
                else if (grepl("sparse", postselection)) {
                  print("Selected covariates are refit with previously estimated weighted ridge prior")
                }
                else {
                  print("Selected covariates are refit with an ordinary ridge prior using newly cross-validated penalty")
                }
            })
        }
        else {
            print("Bondell&Reich marginal credible set used to select covariates")
            ind <- which(penfctr != 0 & beta != 0)
            indpen <- which((beta[penfctr != 0]) != 0)
            if (model == "linear") {
                diagWhalf <- rep(1, n)
            }
            else if (model == "logistic") {
                expminXb <- exp(-X \%*\% c(beta) - intrcpt)
                Pi <- 1/(1 + expminXb)
                diagWhalf <- c(Pi * (1 - Pi))
            }
            else if (model == "cox") {
                expXb <- exp(X \%*\% c(beta))
                h0 <- sapply(1:length(Y[, 1]), function(i) {
                  Y[i, 2]/sum(expXb[Y[, 1] >= Y[i, 1]])
                })
                H0 <- sapply(Y[, 1], function(Ti) {
                  sum(h0[Y[, 1] <= Ti])
                })
                diagWhalf <- c(H0 * expXb)
            }
            Xtilde <- t(t(diagWhalf * X[, ind])/sqrt(lambdap[indpen]))
            svdXtilde <- svd(Xtilde)
            sdBetas <- 1/sqrt(lambdap[indpen]) * sqrt(1 - diag(svdXtilde$v \%*\% 
                (1/(1 + 1/svdXtilde$d^2) * t(svdXtilde$v))))
            s <- sdBetas/min(sdBetas)
            output <- lapply(maxsel2, function(x) {
                fsel <- function(t, maxselec = x) {
                  An <- which(abs(beta[ind]) > t * s)
                  return(length(An) - maxselec)
                }
                maxT <- max(abs(beta[ind])/s) * 2
                t <- uniroot(fsel, interval = c(0, maxT), maxiter = 200, 
                  tol = 10^(-20))$root
                indPost <- ind[which(abs(beta[ind]) > t * s)]
                indPostpen <- which(which(penfctr != 0) \%in\% 
                  indPost)
                indAll <- c(indPost, which(penfctr == 0))
                if (grepl("dense", postselection)) {
                  if (x == maxsel2[1]) 
                    print("Selected covariates are refit with previously estimated weighted ridge prior")
                  lambdaoverall <- sigmahat/tautrgt
                  lam2 <- sigmahat/tautrgt/n * sd_y
                  offset <- rep(0, n)
                  if (length(muhatp) == 1) 
                    muhatp <- rep(muhatp, sum(penfctr != 0))
                  if (!all(muhatp == 0)) 
                    offset <- X[, indAll, drop = F] \%*\% muhatp[indPostpen, 
                      drop = F]
                  Xacc <- X
                  Xacc[, pen] <- as.matrix(X[, pen] \%*\% sparseMatrix(i = 1:length(lambdap), 
                    j = 1:length(lambdap), x = c(1/sqrt(lambdap/lambdaoverall))))
                  if (model == "cox") {
                    glmPost <- glmnet(Xacc[, indAll, drop = F], 
                      Y, alpha = 0, lambda = lam2, family = fml, 
                      offset = offset, standardize = F, penalty.factor = penfctr[indAll], 
                      thresh = 10^-10)
                  }
                  else {
                    glmPost <- glmnet(Xacc[, indAll, drop = F], 
                      Y, alpha = 0, lambda = lam2, family = fml, 
                      offset = offset, intercept = intrcpt, standardize = F, 
                      penalty.factor = penfctr[indAll], thresh = 10^-10)
                  }
                  betaPost <- rep(0, p)
                  betaPost[indAll] <- as.vector(glmPost$beta)
                  betaPost[indPost] <- c(1/sqrt(lambdap[indPostpen]/lambdaoverall)) * 
                    betaPost[indPost] + muhatp[indPostpen]
                }
                else {
                  if (model == "cox") {
                    lambdaGLM <- cv.glmnet(X[, indAll, drop = F], 
                      Y, alpha = 0, family = fml, standardize = F, 
                      penalty.factor = penfctr[indAll])
                    lam2 <- lambdaGLM$lambda.min
                  }
                  else {
                    lambdaGLM <- cv.glmnet(X[, indAll, drop = F], 
                      Y, alpha = 0, family = fml, standardize = F, 
                      penalty.factor = penfctr[indAll])
                    lam2 <- lambdaGLM$lambda.min
                  }
                  if (model == "cox") {
                    glmPost <- glmnet(X[, indAll, drop = F], 
                      Y, alpha = 0, lambda = lam2, family = fml, 
                      offset = offset, standardize = F, penalty.factor = penfctr[indAll], 
                      thresh = 10^-10)
                  }
                  else {
                    glmPost <- glmnet(X[, indAll, drop = F], 
                      Y, alpha = 0, lambda = lam2, family = fml, 
                      intercept = intrcpt, standardize = F, penalty.factor = penfctr[indAll], 
                      thresh = 10^-10)
                  }
                  betaPost <- rep(0, p)
                  betaPost[indAll] <- as.vector(glmPost$beta)
                }
                whichPost <- indPost
                output <- list()
                output$betaPost <- betaPost
                output$whichPost <- whichPost
                output$a0 <- glmPost$a0
                if (!all(is.nan(X2))) {
                  if (model == "linear") {
                    YpredPost <- X2 \%*\% c(betaPost) + glmPost$a0
                    MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                  }
                  if (model == "logistic") {
                    YpredPost <- 1/(1 + exp(-(X2 \%*\% c(betaPost) + 
                      glmPost$a0)))
                    MSEPost <- sum((YpredPost - Y2)^2)/length(Y2)
                    if (any(is.nan(YpredPost))) {
                      browser()
                    }
                  }
                  else if (model == "cox") {
                    expXbPost <- exp(X \%*\% c(betaPost))
                    h0 <- sapply(1:length(Y[, 1]), function(i) {
                      Y[i, 2]/sum(expXbPost[Y[, 1] >= Y[i, 1]])
                    })
                    H0 <- sapply(Y2[, 1], function(Ti) {
                      sum(h0[Y[, 1] <= Ti])
                    })
                    YpredPost <- H0 * exp(X2 \%*\% c(betaPost))
                    MSEPost <- sum((YpredPost - Y2[, 2])^2)/length(Y2[, 
                      2])
                  }
                  output$MSEPost <- MSEPost
                  output$YpredPost <- YpredPost
                }
                return(output)
            })
        }
    }
    res <- list()
    size <- length(maxsel2)
    for (attr in names(output[[1]])) {
        if (attr != "whichPost") {
            res[[attr]] <- matrix(unlist(lapply(output, function(x) {
                x[[attr]]
            })), length(output[[1]][[attr]]), size)
            colnames(res[[attr]]) <- paste("MaxSelec: ", maxsel2, 
                sep = "")
        }
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
