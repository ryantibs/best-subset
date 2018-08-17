n = 100
p = 1000
x = matrix(rnorm(n*p),n,p)
y = rnorm(n)

# Without L1 regularization, cheaper algorithm
t1 = system.time({
  obj1 = L0Learn::L0Learn.fit(x,y,penalty="L0",
                              algorithm="CD",
                              nLambda=50)
})[1]


# With L1 regularization, cheaper algorithm
t2 = system.time({
  obj2 = L0Learn::L0Learn.fit(x,y,penalty="L0L1",
                              algorithm="CD",
                              nGamma=10,nLambda=50)
})[1]

# Without L1 regularization, better algorithm
t3 = system.time({
  obj3 = L0Learn::L0Learn.fit(x,y,penalty="L0",
                              algorithm="CDPSI",
                              nLambda=50)
})[1]

# With L1 regularization, better algorithm
t4 = system.time({
  obj4 = L0Learn::L0Learn.fit(x,y,penalty="L0L1",
                              algorithm="CDPSI",
                              nGamma=10,nLambda=50)
})[1]

cat(t1, "\n", t2, "\n", t3, "\n", t4, "\n")
