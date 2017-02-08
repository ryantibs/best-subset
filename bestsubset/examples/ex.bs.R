# Simulate some simple regression data with the first 5 coefficients
# being nonzero
n = 100
p = 10
xy.obj = sim.xy(n,p,nval=0,ntest=10000,s=5,beta.type=2,snr=0.5)
x = xy.obj$x
y = xy.obj$y

# Run forward stepwise regression for 5 steps
fs.obj = fs(x,y,intercept=FALSE,maxsteps=5,verbose=TRUE)
coef(fs.obj)

# Solve best subset selection for 5 sparsity levels
bs.obj = bs(x,y,intercept=FALSE,k=1:5,verbose=TRUE)
coef(bs.obj)
