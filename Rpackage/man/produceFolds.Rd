\name{produceFolds}
\alias{produceFolds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Produce folds
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Produce folds for cross-validation.
}
\usage{
produceFolds(nsam, outerfold, response, model = c("logistic","cox","other"), 
balance = TRUE, fixedfolds = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsam}{
%%     ~~Describe \code{nsam} here~~
Number of samples
}
  \item{outerfold}{
%%     ~~Describe \code{outerfold} here~~
Number of folds.
}
  \item{response}{
%%     ~~Describe \code{response} here~~
Response data.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
Type of model for the response.
}
  \item{balance}{
%%     ~~Describe \code{balance} here~~
Should folds be balanced in response?
}
  \item{fixedfolds}{
%%     ~~Describe \code{fixedfolds} here~~
Should fixed folds be used?
}
}
