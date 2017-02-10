# Simulate in simple regression setting with the first 5 coefficients
# being nonzero
set.seed(0)
n = 100
p = 20
nval = n
ntest = 10000

# Check for gurobi package
if (!require("gurobi",quietly=TRUE)) {
  stop("Package gurobi not installed (required here)!")
}

# Regression functions: lasso, forward stepwise, and best subset selection
reg.funs = list()
reg.funs[["Lasso"]] = function(x,y) lasso(x,y,intercept=FALSE,nlam=50)
reg.funs[["Stepwise"]] = function(x,y) fs(x,y,intercept=FALSE)
reg.funs[["Best subset"]] = function(x,y) bs(x,y,intercept=FALSE)

# Run the master simulation function, then print results
sim.obj.hisnr = sim.master(n,p,nval,ntest,reg.funs=reg.funs,nrep=10,seed=0,
                     beta.type=2,snr=1,verbose=TRUE)
sim.obj.hisnr

# Repeat, but now for a lower signal-to-noise ratio
sim.obj.losnr = sim.master(n,p,nval,ntest,reg.funs=reg.funs,nrep=10,
                           seed=0,beta.type=2,snr=0.1,verbose=TRUE)
sim.obj.losnr

# Plot simulation results side by side
par(mfrow=c(1,2))
plot(sim.obj.hisnr, main="SNR = 1", legend.pos="topright")
plot(sim.obj.losnr, main="SNR = 0.1", legend.pos="topleft")
