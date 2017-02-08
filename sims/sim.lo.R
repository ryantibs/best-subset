## Low-dimensional simulation, n=500 and p=100
library(bestsubset)

# Set some overall simulation parameters
n = 500; p = 100 # Size of training set, and number of predictors
nval = 500; ntest = 10000 # Size of validation and testing sets
nrep = 50 # Number of repetitions for a given setting
seed = 0 # Random number generator seed
snr.vec = c(seq(0.5,3,length=10), 4:10) # SNRs to consider
rho = 0.8 # Pairwise predictor correlation (for beta.type = 1)
root = "rds/sim.lo"

reg.funs = list()
# Lasso
reg.funs[["Lasso"]] = function(x,y) {
  return(glmnet(x,y,intercept=FALSE))
}
# Forward stepwise
reg.funs[["Stepwise"]] = function(x,y) {
  return(fs(x,y,intercept=FALSE))
}
# Best subset selection
reg.funs[["Best subset"]] = function(x,y) {
  return(bs(x,y,intercept=FALSE))
}


for (beta.type in 1:5) {
  for (snr in snr.vec) {
    file = paste0(root, ".beta", beta.type, ".snr", snr)
    cat("..... NEW SIMULATION .....\n")
    cat("--------------------------\n")
    cat(paste0("File: ", file, "\n\n"))
    
    sim.master(n, p, nval, ntest, reg.funs=reg.funs, nrep=nrep, seed=seed,           
               verbose=TRUE, file=file, rho=rho, beta.type=beta.type, snr=snr)
    
    cat("\n\n")
  }
}
