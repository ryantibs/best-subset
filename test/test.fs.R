library(selectiveInference)
library(bestsubset)

n = 100
p = 20
x = matrix(rnorm(n*p),n,p)
y = rnorm(n)
k = 10

obj1 = selectiveInference::fs(x,y)
obj2 = bestsubset::fs(x,y)

cat(max(abs(obj1$beta-obj2$beta)), "\n")

s = 3.7
cat(max(abs(selectiveInference::coef.fs(obj1,s) -
            bestsubset::coef.fs(obj2,s))), "\n")

x0 = matrix(rnorm(n*p),n,p)
cat(max(abs(selectiveInference::predict.fs(obj1,x0,s) -
            bestsubset::predict.fs(obj2,x0,s))), "\n")
    
