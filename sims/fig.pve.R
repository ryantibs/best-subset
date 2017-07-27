## Degrees of freedom simulation
library(bestsubset)
glmnet.control(fdev=0)

# Set some overall simulation parameters
n = 200; p = 100 # Size of training set, and number of predictors
nval = n # Size of validation set
nrep = 20 # Number of repetitions 
seed = 3 # Random number generator seed
s = 5 # Number of nonzero coefficients
beta.type = 2 # Coefficient type
rho = 0 # Pairwise predictor correlations
snr.vec = seq(0.05,6,length=20) # Signal-to-noise ratios

# Regression function: lasso
reg.funs = list()
reg.funs[["Lasso"]] = function(x,y) lasso(x,y,intercept=FALSE,nlam=50)

sim.list = vector(mode="list",length=length(snr.vec))
for (i in 1:length(snr.vec)) {
  cat("..... NEW SIMULATION .....\n")
  cat("--------------------------\n")
  cat(sprintf("SNR: %0.2f\n\n", snr.vec[i]))
  sim.list[[i]] = sim.master(n, p, nval, reg.funs=reg.funs, nrep=nrep,
                             seed=seed, verbose=TRUE, rho=rho,
                             beta.type=beta.type, snr=snr.vec[i])      
  cat("\n")
}

pve = numeric(length(snr.vec))
for (i in 1:length(snr.vec)) {
  res = tune.and.aggregate(sim.list[[i]], sim.list[[i]]$prop)
  pve[i] = res$z.val.ave
}
snr.vec2 = seq(0,6,length=1000)

dat = data.frame(x=c(snr.vec,snr.vec2),
                 y=c(pve,snr.vec2/(1+snr.vec2)),
                 Type=factor(c(rep("Lasso",length(snr.vec)),
                               rep("Population",length(snr.vec2)))))

ggplot(dat, aes(x=x,y=y,color=Type)) +
  xlab("Signal-to-noise ratio") +
  ylab("Proportion of variance explained") +
  geom_line() + geom_point(aes(shape=Type)) +
  scale_shape_manual(values=c(16,NA)) +
  theme_bw() + theme(legend.just=c(1,0), legend.pos=c(0.95,0.05))
ggsave("fig/pve.pdf", height=3, width=5, device="pdf")
  
