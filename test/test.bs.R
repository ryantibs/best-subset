library(leaps)
library(bestsubset)

n = 100
p = 10
x = matrix(rnorm(n*p),n,p)
y = rnorm(n)

obj1 = leaps::leaps(x,y,nbest=1)
obj2 = bestsubset::bs(x,y,form=1)
obj3 = bestsubset::bs(x,y,form=2)

cat(sum(t(obj1$which) != (obj2$beta[,-1] != 0)), "\n")
cat(sum((obj2$beta != 0) != (obj3$beta != 0)), "\n")
