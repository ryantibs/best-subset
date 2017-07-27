## Degrees of freedom and risk simulations
library(bestsubset)
glmnet.control(fdev=0)

# Set some overall simulation parameters
n = 70; p = 30 # Size of training set, and number of predictors
nval = n # Size of validation set
nrep = 500 # Number of repetitions 
seed = 0 # Random number generator seed
s = 5 # Number of nonzero coefficients
beta.type = 2 # Coefficient type

## Risk simulations
nrep = 20
reg.funs = list()
reg.funs[["Best subset"]] = function(x,y) bs(x,y,intercept=FALSE)
reg.funs[["Forward stepwise"]] = function(x,y) fs(x,y,intercept=FALSE)

# We use a bit of a hack here to get the lasso to return only one solution
# of each size k=0,...,p: for each k, we take the last solution along it
# encounters along the path
reg.funs[["Lasso"]] = function(x,y) {
  out = lasso(x,y,intercept=FALSE,nlam=300)
  class(out) = "lasso2"; return(out)
}
# Now define custom coef and predict functions: where the hacking happens
coef.lasso2 = function(object, s=NULL, gamma=NULL) {
  beta = as.matrix(coef.lasso(object,s,gamma))
  nzs = colSums(beta != 0)
  j = nlam - rev(match(p:0, rev(nzs), NA))
  return(beta[,j])
}
predict.lasso2 = function(object, newx, s=NULL) {
  if (missing(newx)) newx = object$x
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% coef.lasso2(object,s))
}

# Run the master simulation functions at (snr,rho) = (0.35,0.7) (hard setting)
# and (snr,rho) = (0,2) (easy setting)
sim.obj.losnr = sim.master(n,p,nval,reg.funs=reg.funs,nrep=nrep,seed=seed,
                           rho=0.35,s=s,beta.type=beta.type,snr=0.7,verbose=TRUE)
sim.obj.hisnr = sim.master(n,p,nval,reg.funs=reg.funs,nrep=nrep,seed=seed,
                           rho=0.00,s=s,beta.type=beta.type,snr=2.0,verbose=TRUE)

plot(sim.obj.losnr, what="risk", main="SNR=0.7, Cor=0.35", legend=FALSE,
     make.pdf=TRUE, fig.dir="fig", file.name="snr.lo", h=4, w=4)
plot(sim.obj.hisnr, what="risk", main="SNR=2.0, Cor=0.00", make.pdf=TRUE,
     fig.dir="fig", file.name="snr.hi", h=4, w=5.5)

## Degrees of freedom simulation
# Note this is a "hand-made" simulation (rather than using built-in functions
# like above) because we want to keep x fixed for the df calculations
set.seed(seed)
xy.obj = sim.xy(n,p,nval,rho=0.35,s=s,beta.type=beta.type,snr=0.7)
x = xy.obj$x
y = xy.obj$y
mu = as.numeric(x %*% xy.obj$beta)
sigma = xy.obj$sigma
nlam = 300
nrel = 9

ip.las = matrix(0,nrep,(p+1)*nrel)
ip.fs = ip.bs = matrix(0,nrep,p+1)

for (r in 1:nrep) {
  cat(r,"... ")
  eps = rnorm(n)*sigma
  y = mu + eps
  beta.las = coef(lasso(x,y,intercept=FALSE,nlam=nlam,nrel=nrel))
  nzs.las = colSums(beta.las != 0)[0:(nlam-1)*nrel+1]
  j = nlam - rev(match(p:0, rev(nzs.las), NA))
  ind = rep(j*nrel, each=nrel) + rep(1:nrel, p+1)
  yhat.las = (x %*% beta.las)[,ind]
  ip.las[r,] = colSums(yhat.las * eps)
  
  yhat.fs = predict(fs(x,y,intercept=FALSE))
  yhat.bs = predict(bs(x,y,intercept=FALSE))
  ip.fs[r,] = colSums(yhat.fs * eps)
  ip.bs[r,] = colSums(yhat.bs * eps)
}

df.las = colMeans(ip.las, na.rm=TRUE) / sigma^2
df.las = matrix(df.las, p+1, nrel, byrow=TRUE)
df.fs = colMeans(ip.fs, na.rm=TRUE) / sigma^2
df.bs = colMeans(ip.bs, na.rm=TRUE) / sigma^2

save(list=ls(),file="sim.df.rda")

##############################
# Run the code below to reproduce the df figure without rerunning the sims
library(bestsubset)
load("rds/sim.df.rda")

# Plot the results
dat = data.frame(x=rep(0:p,3),
                 y=c(df.bs,df.fs,df.las[,1]),
                 Method=factor(rep(c("Best subset","Forward stepwise","Lasso"),
                                    rep(p+1,3))))

ggplot(dat, aes(x=x,y=y,color=Method)) +
  xlab("Number of nonzero coefficients") +
  ylab("Degrees of freedom") +
  geom_line(lwd=0.5, color="black", linetype=3, aes(x,x)) +
  geom_line(lwd=1) + geom_point(pch=19) +
  theme_bw() + theme(legend.just=c(1,0), legend.pos=c(0.95,0.05)) 
ggsave("fig/df1.pdf", height=4, width=4, device="pdf")

dat = data.frame(x=rep(0:p,5),
                 y=c(df.bs,df.fs,df.las[,1],df.las[,5],df.las[,9]),
                 Method=factor(rep(c("Best subset","Forward stepwise","Lasso",
                                      "Relaxed lasso: 0.5","Relaxed lasso: 0"),
                                    rep(p+1,5))))

ggplot(dat, aes(x=x,y=y,color=Method)) +
  xlab("Number of nonzero coefficients") +
  ylab("Degrees of freedom") +
  geom_line(lwd=0.5, color="black", linetype=3, aes(x,x)) +
  geom_line(lwd=1) + geom_point(pch=19) + 
  theme_bw() + theme(legend.just=c(1,0), legend.pos=c(0.95,0.05)) 
ggsave("fig/df2.pdf", height=4, width=4, device="pdf")

