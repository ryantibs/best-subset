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

# Plot simulation results, excluding relaxed lasso 
par(mfrow=c(1,2))
plot(sim.obj.hisnr, method.nums=1:3, main="SNR = 1", legend.pos="topright")
plot(sim.obj.losnr, method.nums=1:3, main="SNR = 0.1", legend.pos="topleft")

# Plot simulation results, including relaxed lasso (it looks a bit crazy)
par(mfrow=c(1,2))
plot(sim.obj.hisnr, main="SNR = 1", legend.pos="topright")
plot(sim.obj.losnr, main="SNR = 0.1", legend.pos="topleft")
