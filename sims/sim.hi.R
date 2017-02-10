## High-dimensional simulation, n=50 and p=1000
library(bestsubset)
library(glmnet)

# Set some overall simulation parameters
n = 50; p = 1000 # Size of training set, and number of predictors
nval = n; ntest = 10000 # Size of validation and testing sets
nrep = 10 # Number of repetitions for a given setting
seed = 0 # Random number generator seed
type.vec = 1:5 # Simulation settings to consider
snr.vec = c(seq(0.01,1,length=10), 2:4) # SNRs to consider
rho = 0.8 # Pairwise predictor correlation (for beta.type = 1)
stem = "sim.lo"

# Regression functions: lasso, forward stepwise, and best subset selection
reg.funs = list()
reg.funs[["Lasso"]] = function(x,y) lasso(x,y,intercept=FALSE,nlam=50)
reg.funs[["Stepwise"]] = function(x,y) fs(x,y,intercept=FALSE)
reg.funs[["Best subset"]] = function(x,y) bs(x,y,intercept=FALSE)

file.list = c()
file.name = c()
for (beta.type in type.vec) {
  for (snr in snr.vec) {
    name = paste0(stem, ".beta", beta.type)
    file = paste0("rds/", name, ".snr", snr, ".rds")
    cat("..... NEW SIMULATION .....\n")
    cat("--------------------------\n")
    cat(paste0("File: ", file, "\n\n"))
    
    sim.master(n, p, nval, ntest, reg.funs=reg.funs, nrep=nrep, seed=seed,           
               verbose=TRUE, file=file, rho=rho, beta.type=beta.type, snr=snr)

    file.list = c(file.list, file)
    file.name = c(file.name, name)
    cat("\n")
  }
}

plot.many.sims(file.list, grouping=rep(type.vec, each=length(snr.vec)),
               snr.vec=snr.vec, fig.dir="fig", file.name=file.name)

