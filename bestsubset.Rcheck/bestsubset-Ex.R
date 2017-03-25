pkgname <- "bestsubset"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bestsubset')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("bs")
### * bs

flush(stderr()); flush(stdout())

### Name: bs
### Title: Best subset selection.
### Aliases: bs

### ** Examples

# Simulate some simple regression data with the first 5 coefficients
# being nonzero
set.seed(3)
n = 100
p = 20
ntest = 10000
xy.obj = sim.xy(n,p,nval=0,ntest=ntest,s=5,beta.type=2,snr=1)
x = xy.obj$x
y = xy.obj$y
xtest = xy.obj$xtest
mutest = xy.obj$mutest

# Run forward stepwise regression for 8 steps
fs.obj = fs(x,y,intercept=FALSE,maxsteps=8,verbose=TRUE)
fs.beta = coef(fs.obj)
fs.supp = apply(fs.beta != 0, 2, which)

# Solve best subset selection for 8 sparsity levels
bs.obj = bs(x,y,intercept=FALSE,k=1:8,verbose=TRUE)
bs.beta = coef(bs.obj)
bs.supp = apply(bs.beta != 0, 2, which)

# Compare supports of the solutions with 5 and 8 variables
fs.supp[[5]]; bs.supp[[5]]
fs.supp[[8]]; bs.supp[[8]]

# Predict on test data and record risk
fs.pred = predict(fs.obj,newx=xtest)
bs.pred = predict(fs.obj,newx=xtest)
colMeans((fs.pred - mutest)^2)
colMeans((bs.pred - mutest)^2)



cleanEx()
nameEx("fs")
### * fs

flush(stderr()); flush(stdout())

### Name: fs
### Title: Forward stepwise regression.
### Aliases: fs

### ** Examples

# Simulate some simple regression data with the first 5 coefficients
# being nonzero
set.seed(3)
n = 100
p = 20
ntest = 10000
xy.obj = sim.xy(n,p,nval=0,ntest=ntest,s=5,beta.type=2,snr=1)
x = xy.obj$x
y = xy.obj$y
xtest = xy.obj$xtest
mutest = xy.obj$mutest

# Run forward stepwise regression for 8 steps
fs.obj = fs(x,y,intercept=FALSE,maxsteps=8,verbose=TRUE)
fs.beta = coef(fs.obj)
fs.supp = apply(fs.beta != 0, 2, which)

# Solve best subset selection for 8 sparsity levels
bs.obj = bs(x,y,intercept=FALSE,k=1:8,verbose=TRUE)
bs.beta = coef(bs.obj)
bs.supp = apply(bs.beta != 0, 2, which)

# Compare supports of the solutions with 5 and 8 variables
fs.supp[[5]]; bs.supp[[5]]
fs.supp[[8]]; bs.supp[[8]]

# Predict on test data and record risk
fs.pred = predict(fs.obj,newx=xtest)
bs.pred = predict(fs.obj,newx=xtest)
colMeans((fs.pred - mutest)^2)
colMeans((bs.pred - mutest)^2)



cleanEx()
nameEx("sim.master")
### * sim.master

flush(stderr()); flush(stdout())

### Name: sim.master
### Title: Master function for running simulations.
### Aliases: sim.master

### ** Examples

# Simulate in simple regression setting with the first 5 coefficients
# being nonzero
set.seed(0)
n = 100
p = 20
nval = n

# Check for gurobi package
if (!require("gurobi",quietly=TRUE)) {
  stop("Package gurobi not installed (required here)!")
}

# Regression functions: lasso, forward stepwise, and best subset selection
reg.funs = list()
reg.funs[["Lasso"]] = function(x,y) lasso(x,y,intercept=FALSE,nlam=50)
reg.funs[["Stepwise"]] = function(x,y) fs(x,y,intercept=FALSE)
reg.funs[["Best subset"]] = function(x,y) bs(x,y,intercept=FALSE)
reg.funs[["Relaxed lasso"]] = function(x,y) lasso(x,y,intercept=FALSE,
                                                  nrelax=5,nlam=50)

# Run the master simulation function, for two different SNRs 
sim.obj.hisnr = sim.master(n,p,nval,reg.funs=reg.funs,nrep=10,seed=0,
                           beta.type=2,s=5,snr=1,verbose=TRUE)
sim.obj.losnr = sim.master(n,p,nval,reg.funs=reg.funs,nrep=10,seed=0,
                           beta.type=2,s=5,snr=0.1,verbose=TRUE)

# Print simulation results
sim.obj.hisnr
sim.obj.losnr

# Plot simulation results, excluding relaxed lasso (it looks a bit crazy)
par(mfrow=c(1,2))
plot(sim.obj.hisnr, method.nums=1:3, main="SNR = 1", legend.pos="topright")
plot(sim.obj.losnr, method.nums=1:3, main="SNR = 0.1", legend.pos="topleft")

# Plot simulation results, including relaxed lasso (it looks a bit crazy)
par(mfrow=c(1,2))
plot(sim.obj.hisnr, main="SNR = 1", legend.pos="topright")
plot(sim.obj.losnr, main="SNR = 0.1", legend.pos="topleft")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("sim.xy")
### * sim.xy

flush(stderr()); flush(stdout())

### Name: sim.xy
### Title: Predictors and responses generation.
### Aliases: sim.xy

### ** Examples

# Simulate some simple regression data with the first 5 coefficients
# being nonzero
set.seed(3)
n = 100
p = 20
ntest = 10000
xy.obj = sim.xy(n,p,nval=0,ntest=ntest,s=5,beta.type=2,snr=1)
x = xy.obj$x
y = xy.obj$y
xtest = xy.obj$xtest
mutest = xy.obj$mutest

# Run forward stepwise regression for 8 steps
fs.obj = fs(x,y,intercept=FALSE,maxsteps=8,verbose=TRUE)
fs.beta = coef(fs.obj)
fs.supp = apply(fs.beta != 0, 2, which)

# Solve best subset selection for 8 sparsity levels
bs.obj = bs(x,y,intercept=FALSE,k=1:8,verbose=TRUE)
bs.beta = coef(bs.obj)
bs.supp = apply(bs.beta != 0, 2, which)

# Compare supports of the solutions with 5 and 8 variables
fs.supp[[5]]; bs.supp[[5]]
fs.supp[[8]]; bs.supp[[8]]

# Predict on test data and record risk
fs.pred = predict(fs.obj,newx=xtest)
bs.pred = predict(fs.obj,newx=xtest)
colMeans((fs.pred - mutest)^2)
colMeans((bs.pred - mutest)^2)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
