## plots for high-dimensional simulation, s=10
library(bestsubset)
n = 100; p = 1000
file.list = system(paste0("ls ","rds/hi10/*.rds"),intern=TRUE)
method.nums = c(7,2,1,3)
method.names = c("Best subset","Forward stepwise","Lasso","Relaxed lasso")


plot.from.file(file.list, what="error", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="fig/hi10",
               file.name="hi10.err")

plot.from.file(file.list, what="prop", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="fig/hi10",
               file.name="hi10.prop")

plot.from.file(file.list, what="nonzero", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="fig/hi10",
               file.name="hi10.nzs")

plot.from.file(file.list, what="F", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="fig/hi10",
               file.name="hi10.F")
