library(selectiveInference)
library(bestsubset)

n = 100
p = 20
x = matrix(rnorm(n*p),n,p)
y = rnorm(n)

###

obj1 = selectiveInference::fs(x,y,intercept=FALSE)
obj2 = bestsubset::fs(x,y,intercept=FALSE)

s = 3.7
cat(max(abs(selectiveInference::coef.fs(obj1,s+1) -
            bestsubset::coef.fs(obj2,s))), "\n")

x0 = matrix(rnorm(n*p),n,p)
cat(max(abs(selectiveInference::predict.fs(obj1,x0,s+1) -
            bestsubset::predict.fs(obj2,x0,s))), "\n")
    
cat(max(abs(selectiveInference::coef.fs(obj1,s=p+1) - 
            bestsubset::coef.fs(obj2,s=p))))

###

obj1 = selectiveInference::fs(x,y,intercept=TRUE)
obj2 = bestsubset::fs(x,y,intercept=TRUE)

s = 3.7
cat(max(abs(selectiveInference::coef.fs(obj1,s+1) -
            bestsubset::coef.fs(obj2,s))), "\n")

x0 = matrix(rnorm(n*p),n,p)
cat(max(abs(selectiveInference::predict.fs(obj1,x0,s+1) -
            bestsubset::predict.fs(obj2,x0,s))), "\n")
    
cat(max(abs(selectiveInference::coef.fs(obj1,s=p+1) - 
            bestsubset::coef.fs(obj2,s=p))))
