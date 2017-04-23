library(glmnet)
library(bestsubset)

n = 100
p = 20
x = matrix(rnorm(n*p),n,p)
y = rnorm(n)

###

obj1 = glmnet::glmnet(x,y,intercept=FALSE,nlambda=10)
obj2 = bestsubset::lasso(x,y,intercept=FALSE,lambda=obj1$lambda)

cat(max(abs(glmnet::coef.glmnet(obj1)[-1,] -
            bestsubset::coef.lasso(obj2))), "\n")

x0 = matrix(rnorm(n*p),n,p)
cat(max(abs(glmnet::predict.glmnet(obj1,x0) -
            bestsubset::predict.lasso(obj2,x0))), "\n")

###

obj1 = glmnet::glmnet(x,y,intercept=TRUE,nlambda=10)
obj2 = bestsubset::lasso(x,y,intercept=TRUE,lambda=obj1$lambda)

cat(max(abs(glmnet::coef.glmnet(obj1) -
            bestsubset::coef.lasso(obj2))), "\n")

x0 = matrix(rnorm(n*p),n,p)
cat(max(abs(glmnet::predict.glmnet(obj1,x0) -
            bestsubset::predict.lasso(obj2,x0))), "\n")
