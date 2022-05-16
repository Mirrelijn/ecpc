# ecpc
Use the ecpc method to learn from co-data.
The R-package may be downloaded from CRAN as follows:

install.packages("ecpc")

or from this github-page;

library(devtools); install_github("Mirrelijn/ecpc/Rpackage")

A [cheat sheet]("Cheatsheet.pdf") for the functions is:
![Cheatsheet](https://github.com/Mirrelijn/ecpc/blob/master/Cheatsheet_ecpc.png?raw=true)

-------

The file "Demo.R" contains a short demo of ecpc on toy data.

The folders "Simulation study", "Application miRNA data" and "Application methylation data" contain the code and data for reproduction of the analyses in the paper.

The R-package ecpc is short for: Empirical bayes Co-data learnt Prediction and Covariate selection. 
The method learns from so-called co-data to improve prediction and covariate selection in high-dimensional data. Co-data complements the main data by grouping covariates into groups that may differ in predictive strength. Each group obtains a group weight to quantify this predictive strength, or average effect size. The group weights are estimated using an empirical Bayes method of moments with an extra later of hypershrinkage. Different types of hypershrinkage may be used suitable for the type of co-data, e.g. to counter overfitting for many groups or non-informative co-data and to incorporate structure on group level.

Currently, the software allows for:

-Linear, logistic and Cox survival response;

-Multiple co-data, possibly varying in type: non-overlapping groups, overlapping groups, hierarchical groups and continuous co-data;

-Different types of hypershrinkage: ridge for correlated or many groups, lasso for group selection, hierarchical lasso for selection of hierarchical groups and generalised ridge shrinkage for generalised additive co-data models or shape-constrained additive co-data models;

-Inclusion of unpenalised covariates (e.g. clinical covariates like age and gender);

-Posterior selection for a user-defined number of variables.
