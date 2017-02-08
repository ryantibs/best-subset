# Best Subset Selection Simulations
## Trevor Hastie, Rob Tibshirani, Ryan Tibshirani
### Maintained by Ryan Tibshirani

This project contains tools for running simulations that compare best subset
selection in regression to other common sparse regression estimators such as
the lasso and forward stepwise. 

The simulation setup is based on the paper: Bertsimas, King, Mazumder (2016), 
"Best subset selection via a modern optimization lens", as is the mixed integer
programming formulation used for solving best subset selection.

### Install the R package

To install the conformalInference R package directly from github, run the
following in R:

```{r}
library(devtools)
install_github(repo="ryantibs/best-subset", subdir="bestsubset")
```  