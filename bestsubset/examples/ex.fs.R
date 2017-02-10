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
