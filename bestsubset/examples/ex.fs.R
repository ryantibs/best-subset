# Simulate some simple regression data with the first 5 coefficients
# being nonzero
set.seed(0)
n = 100
p = 20
xy.obj = sim.xy(n,p,nval=0,s=5,beta.type=2,snr=1)
x = xy.obj$x
y = xy.obj$y

# Run forward stepwise regression for 8 steps
fs.obj = fs(x,y,intercept=FALSE,maxsteps=8,verbose=TRUE)
fs.beta = coef(fs.obj)
fs.supp = apply(fs.beta != 0, 2, which)

# Solve best subset selection for 8 sparsity levels
bs.obj = bs(x,y,intercept=FALSE,k=0:8,verbose=TRUE)
bs.beta = coef(bs.obj)
bs.supp = apply(bs.beta != 0, 2, which)

# Compare supports of the solutions with 5 and 8 variables
fs.supp[[6]]; bs.supp[[6]]
fs.supp[[9]]; bs.supp[[9]]

# Predict on test data and record test error
ntest = 10000
xy.obj.test = sim.xy(ntest,p,nval=0,s=5,beta.type=2,snr=1)
xtest = xy.obj.test$x
ytest = xy.obj.test$y

fs.pred = predict(fs.obj,newx=xtest)
bs.pred = predict(bs.obj,newx=xtest)
colMeans((fs.pred - ytest)^2)
colMeans((bs.pred - ytest)^2)
