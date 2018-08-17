n = 100
p = 10
x = matrix(rnorm(n*p),n,p)
y = rnorm(n)

obj1 = leaps::leaps(x,y,nbest=1)
obj2 = bestsubset::bs(x,y,form=1)
obj3 = bestsubset::bs(x,y,form=2)
obj4 = L0Learn::L0Learn.fit(x,y)

cat(sum(t(obj1$which) != (obj2$beta[,-1] != 0)), "\n")
cat(sum((obj2$beta != 0) != (obj3$beta != 0)), "\n")

cat(all(unlist(obj4$converged)), "\n")
k = obj4$suppSize[[1]][median(1:p)]
cat(sum(obj1$which[k,] !=
        (coef(obj4, lambda=obj4$lambda[[1]][k]) != 0)[-1]), "\n")
