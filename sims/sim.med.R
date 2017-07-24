## Medium-dimensional simulation, n=500 and p=100
library(bestsubset)
glmnet.control(fdev=0)

# Set some overall simulation parameters
n = 500; p = 100 # Size of training set, and number of predictors
nval = n # Size of validation set
nrep = 10 # Number of repetitions for a given setting
seed = 0 # Random number generator seed
type.vec = c(1:3,5) # Simulation settings to consider
rho.vec = c(0,0.35,0.7) # Pairwise predictor correlations
snr.vec = exp(seq(log(0.05),log(6),length=10)) # Signal-to-noise ratios 
stem = paste0("sim.n",n,".p",p)

# Regression functions: lasso, forward stepwise, and best subset selection
reg.funs = list()
reg.funs[["Lasso"]] = function(x,y) lasso(x,y,intercept=FALSE,nlam=100)
reg.funs[["Forward stepwise"]] = function(x,y) fs(x,y,intercept=FALSE,
                                                  max=50)
reg.funs[["Best subset"]] = function(x,y) bs(x,y,intercept=FALSE,k=0:50,
                                             time.limit=180)
reg.funs[["Relaxed lasso"]] = function(x,y) lasso(x,y,intercept=FALSE,
                                                  nrelax=10,nlam=100)

file.list = c() # Vector of files for the saved rds files
for (beta.type in type.vec) {
  for (rho in rho.vec) {
    name = paste0(stem, ".beta", beta.type, sprintf(".rho%0.2f", rho))
    for (snr in snr.vec) {
      file = paste0("rds/", name, ".snr", round(snr,2), ".rds")
      cat("..... NEW SIMULATION .....\n")
      cat("--------------------------\n")
      cat(paste0("File: ", file, "\n\n"))
      
      sim.master(n, p, nval, reg.funs=reg.funs, nrep=nrep, seed=seed,
                 verbose=TRUE, file=file, rho=rho, beta.type=beta.type, snr=snr)
      
      file.list = c(file.list, file)
      cat("\n")
    }
  }
}

##############################
# Run the code below to reproduce the figures without rerunning the sims

library(bestsubset)
n = 500; p = 100
file.list = system(paste0("ls rds/sim.n",n,".p",p,".*.rds"),intern=TRUE)
method.nums = c(3,2,1,4)
method.names = c("Best subset","Forward stepwise","Lasso","Relaxed lasso")

# Validation tuning
plot.from.file(file.list, what="risk", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",5), make.pdf=TRUE,
               fig.dir="fig/val",
               file.name=paste0("sim.n",n,".p",p,".val.risk.rel"))

plot.from.file(file.list, what="error", rel.to=3, tuning="val",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",5), make.pdf=TRUE,
               fig.dir="fig/val", 
               file.name=paste0("sim.n",n,".p",p,".val.err.rel"))

plot.from.file(file.list, what="prop", tuning="val",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",5), make.pdf=TRUE,
               fig.dir="fig/val", 
               file.name=paste0("sim.n",n,".p",p,".val.prop"))

plot.from.file(file.list, what="nonzero", tuning="val",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",5), make.pdf=TRUE,
               fig.dir="fig/val", 
               file.name=paste0("sim.n",n,".p",p,".val.nzs"))

# Oracle tuning
plot.from.file(file.list, what="risk", rel.to=NULL, tuning="ora",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",5), make.pdf=TRUE,
               fig.dir="fig/ora", 
               file.name=paste0("sim.n",n,".p",p,".ora.risk.rel"))

plot.from.file(file.list, what="error", rel.to=3, tuning="ora",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",5), make.pdf=TRUE,
               fig.dir="fig/ora", 
               file.name=paste0("sim.n",n,".p",p,".ora.err.rel"))

plot.from.file(file.list, what="prop", tuning="ora",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",5), make.pdf=TRUE,
               fig.dir="fig/ora", 
               file.name=paste0("sim.n",n,".p",p,".ora.prop"))

plot.from.file(file.list, what="nonzero", tuning="ora",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",5), make.pdf=TRUE,
               fig.dir="fig/ora", 
               file.name=paste0("sim.n",n,".p",p,".ora.nzs"))
