\name{ecpc}
\alias{ecpc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fit adaptive multi-group ridge GLM with hypershrinkage
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Fits a generalised linear (linear, logistic) or Cox survival model, penalised with adaptive multi-group ridge penalties.
The multi-group penalties correspond to groups of covariates defined by (multiple) co-data sources. 
Group hyperparameters are estimated with an empirical Bayes method of moments, penalised with an extra level of hypershrinkage.
Various types of hypershrinkage may be used for various co-data, including overlapping groups, hierarchical groups and continuous co-data.
}
\usage{
ecpc(Y, X, groupsets, groupsets.grouplvl = NULL, hypershrinkage, 
unpen = NULL, intrcpt = TRUE, model=c("linear","logistic","cox"), 
postselection = "elnet,dense", maxsel = 10, 
lambda = NULL, fold = 10, sigmasq = NaN, w = NaN, 
nsplits = 100, weights = TRUE, profplotRSS = FALSE, Y2 = NaN, X2 = NaN, 
compare = TRUE, mu = FALSE, normalise = FALSE, silent = FALSE,
datablocks = NULL)%,standardise_Y=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
Response data; n-dimensional vector (n: number of samples) for linear and logistic outcomes, or \code{\link{Surv}} object for Cox survival.
}
  \item{X}{
%%     ~~Describe \code{X} here~~
Observed data; (nxp)-dimensional matrix (p: number of covariates) with each row the observed high-dimensional feature vector of a sample.
}
  \item{groupsets}{
%%     ~~Describe \code{groupsets} here~~
Co-data group sets; list with m (m: number of group sets) group sets. Each group set is a list of all groups in that set. Each group is a vector containing the indices of the covariates in that group.
}
  \item{groupsets.grouplvl}{
%%     ~~Describe \code{groupsets.grouplvl} here~~
(optional) Group sets on group level used in hypershrinkage; list of m elements (corresponding to 'groupsets'), with NULL if there is no structure on group level, or with a list of groups containing the indices of groups of covariates in that group. May be used for hierarchical groups and to adaptively discretise continuous co-data, see \code{\link{obtainHierarchy}}.
}
  \item{hypershrinkage}{
%%     ~~Describe \code{hypershrinkage} here~~
Type of shrinkage that is used on the group level; vector of m strings indicating the shrinkage type (or penalty) that is used for each of the m group sets. String may be of the simple form "type1", or "type1,type2", in which type1 is used to select groups and type2 to estimate the group weights of the selected groups. Possible hypershrinkage types are:

c("none","ridge","lasso","hierLasso","lasso,ridge","hierLasso,ridge"); 

"none" for no hypershrinkage, "ridge" (default), "lasso" and "hierLasso" (hierarchical lasso using a latent overlapping group lasso penalty) for group selection possibly be combined with ridge shrinkage.
}
  \item{unpen}{
%%     ~~Describe \code{unpen} here~~
Unpenalised covariates; vector with indices of covariates that should not be penalised.
}
  \item{intrcpt}{
%%     ~~Describe \code{intrcpt} here~~
Should an intercept be included? Included by default for linear and logistic, excluded for Cox for which the baseline hazard is estimated.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
Type of model for the response; linear, logistic or Cox.
}
  \item{postselection}{
%%     ~~Describe \code{postselection} here~~
Type of posterior selection method used to obtain a parsimonious model of maxsel covariates, or FALSE if no parsimonious model is needed. Possible options are "elnet,dense" (default), "elnet,sparse", "BRmarginal,dense", "BRmarginal,sparse" or "DSS".
}
  \item{maxsel}{
%%     ~~Describe \code{maxsel} here~~
Maximum number of covariates to be selected a posteriori, in addition to all unpenalised covariates. If maxsel is a vector, multiple parsimonious models are returned.
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
Global ridge penalty; if given, numeric value to fix the global ridge penalty and equivalently, the global prior variance. When not given, for linear, by default "ML" is used for estimation for maximum marginal likelihood estimation and "CV" for other models for cross-validation.
}
  \item{fold}{
%%     ~~Describe \code{fold} here~~
Number of folds used in inner cross-validation to estimate global ridge penalty lambda.
}
  \item{sigmasq}{
%%     ~~Describe \code{sigmasq} here~~
(linear model only) If given, noise level is fixed (Y~N(X*beta,sd=sqrt(sigmasq))).
}
  \item{w}{
%%     ~~Describe \code{w} here~~
Group set weights: m-dimensional vector. If given, group set weights are fixed.
}
  \item{nsplits}{
%%     ~~Describe \code{nsplits} here~~
Number of splits used in the Residual Sum of Squares (RSS) criterion to estimate the optimal hyperlambda.
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
Should weights be used in hypershrinkage to correct for group size (default TRUE)?
}
  \item{profplotRSS}{
  Should a profile plot of the residual sum of squares (RSS) criterium be shown?
  }
  \item{Y2}{
%%     ~~Describe \code{Y2} here~~
(optional) Independent response data to compare with predicted response.
}
  \item{X2}{
%%     ~~Describe \code{X2} here~~
(optional) Independent observed data for which response is predicted.
}
  \item{compare}{
%%     ~~Describe \code{compare} here~~
Should an ordinary ridge model be fitted to compare with?
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
Should group prior means be included (default FALSE)?
}
  \item{normalise}{
%%     ~~Describe \code{normalise} here~~
Should group variances be normalised to sum to 1 (default FALSE)?
}
  \item{silent}{
%%     ~~Describe \code{normalise} here~~
Should output messages be suppressed (default FALSE)?
}
  \item{datablocks}{
%%     ~~Describe \code{normalise} here~~
(optional) for multiple data types, the corresponding blocks of data may be given in datablocks; a list of B vectors of the indices of covariates in `X' that belong to each of the B data blocks. Unpenalised covariates should not be given as seperate block, but can be omitted or included in blocks with penalised covariates. Each datatype obtains a datatype-specific `tauglobal' as in multiridge.
}
%  \item{standardise_Y}{
%%     ~~Describe \code{normalise} here~~
%Should linear response be standardised (default FALSE)? Regression coefficients are scaled back before returning the %output.
%}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list with the following elements:
  \item{beta}{Estimated regression coefficients; p-dimensional vector.}
  \item{intercept}{If included, the estimated intercept; scalar.}
  \item{tauglobal}{Estimated global prior variance; scalar (or vector with datatype-specific global prior variances when multiple `datablocks' are given).)}
  \item{gammatilde}{Estimated group weights before truncating negative weights to 0; vector of dimension the total number of groups.}
  \item{gamma}{Final estimated group weights; vector of dimension the total number of groups.}
  \item{w}{Estimated group set weights; m-dimensional vector.}
  \item{penalties}{Estimated multi-group ridge penalties; p-dimensinional vector.}
  \item{hyperlambdas}{Estimated hyperpenalty parameters used in hypershrinkage; m-dimensional vector.}
  \item{Ypred}{If independent test set 'X2' is given, predictions for the test set.}
  \item{MSEecpc}{If independent test set 'X2', 'Y2' is given, mean squared error of the predictions.}
  \item{sigmahat}{(linear model) Estimated sigma^2.}
If 'compare'=TRUE, ordinary ridge estimates and predictions are given. If in addition multiple `datablocks' are given, the estimates and predictions for multiridge penalty are given;
  \item{betaridge}{Estimated regression coefficients for ordinary ridge (or multiridge) penalty.}
  \item{interceptridge}{Estimated intercept for ordinary ridge (or multiridge) penalty.}
  \item{lambdaridge}{Estimated (multi)ridge penalty.}
  \item{Ypredridge}{If independent test set 'X2' is given, ordinary ridge (or multiridge) predictions for the test set.}
  \item{MSEridge}{If independent test set 'X2', 'Y2' is given, mean squared error of the ordinary ridge (or multiridge) predictions.}
If posterior selection is performed;
  \item{betaPost}{Estimated regression coefficients for parsimonious models. If 'maxsel' is a vector, 'betaPost' is a matrix with each column the vector estimate corresponding to the maximum number of selected covariates given in 'maxsel'.}
  \item{interceptPost}{Estimated intercept coefficient for parsimonious models.}
  \item{YpredPost}{If independent test set 'X2' is given, posterior selection model predictions for the test set.}
  \item{MSEPost}{If independent test set 'X2', 'Y2' is given, mean squared error of the posterior selection model predictions.}
}
\references{
%% ~put references to the literature/web site here ~
Van Nee, Mirrelijn M., Lodewyk FA Wessels, and Mark A. van de Wiel. "Flexible co-data learning for high-dimensional prediction." arXiv preprint arXiv:2005.04010 (2020).

Van de Wiel, Mark A., Mirrelijn M. van Nee, and Armin Rauschenberger. "Fast cross-validation for multi-penalty ridge regression." arXiv preprint arXiv:2005.09301 (2020).
}
\author{
%%  ~~who you are~~
Mirrelijn van Nee, Lodewyk Wessels, Mark van de Wiel
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\donttest{ 
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

###########################
# Make co-data group sets #
###########################
#Group set 1: G random groups
G <- 5 #number of groups
#sample random categorical co-data:
categoricalRandom <- as.factor(sample(1:G,p,TRUE))
#make group set, i.e. list with G groups:
groupsetRandom <- createGroupset(categoricalRandom)

#Group set 2: informative hierarchical group set
continuousCodata <- abs(Dat$beta) #use the magnitude of beta as continuous co-data
#Use adaptive discretisation to find a good discretisation of the continuous co-data;
# discretise in groups of covariates of various sizes:
groupsetHierarchical <- splitMedian(values=continuousCodata,index = 1:p,
                        minGroupSize = 50,split="both") 
# and obtain group set on group level that defines the hierarchy:
hierarchy.grouplevel <- obtainHierarchy(groupset = groupsetHierarchical) 
#visualise hierarchical groups:
#visualiseGroupset(Groupset = groupsetHierarchical,groupset.grouplvl = hierarchy.grouplevel) 
 
############ 
# Fit ecpc #
############

#fit ecpc for the two group sets, with ridge hypershrinkage for group set 1, 
# and hierarchical lasso and ridge for group set 2.
tic<-proc.time()[[3]]
fit <- ecpc(Y=Dat$Y,X=Dat$Xctd,groupsets=list(groupsetRandom,groupsetHierarchical),
           groupsets.grouplvl=list(NULL,hierarchy.grouplevel),
           hypershrinkage=c("ridge","hierLasso,ridge"),
           model="linear",maxsel=c(5,10,15,20),
           Y2=Dat$Y2,X2=Dat$X2ctd)
toc <- proc.time()[[3]]-tic

fit$tauglobal #estimated global prior variance
fit$gamma #estimated group weights (concatenated for the group sets)
fit$w #estimated group set weights
summary(fit$beta) #estimated regression coefficients
summary(fit$betaPost) #estimated regression coefficients after posterior selection

c(fit$MSEecpc,fit$MSEridge) #mean squared error on test set for ecpc and ordinary ridge
fit$MSEPost #MSE on the test set of ecpc after posterior selection


###################################
# Fit ecpc for multiple datatypes #
###################################
rankBeta<-order(abs(Dat$beta)) #betas ranked in order of magnitude
 
#with multiple datatypes (given in datablocks) and informative groups
fit2 <- ecpc(Y=Dat$Y,X=Dat$Xctd[,rankBeta],groupsets=list(list(1:75,76:150,151:225,226:300)),
            groupsets.grouplvl=list(NULL),
            hypershrinkage=c("none"),
            model="linear",maxsel=c(5,10,15,20),
            Y2=Dat$Y2,X2=Dat$X2ctd[,rankBeta],
            datablocks = list(1:floor(p/2),(floor(p/2)+1):p))
}

}
