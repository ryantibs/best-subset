## plots for high-dimensional simulation, s=5
library(bestsubset)
n = 500; p = 100
file.list = system(paste0("ls ","rds/med/*.rds"),intern=TRUE)
method.nums = c(7,2,1,3)
method.names = c("Best subset","Forward stepwise","Lasso","Relaxed lasso")


plot.from.file(file.list, what="error", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="fig/med",
               file.name="med.err")

plot.from.file(file.list, what="prop", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="fig/med",
               file.name="med.prop")

plot.from.file(file.list, what="nonzero", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="fig/med",
               file.name="med.nzs")g

plot.from.file(file.list, what="F", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="fig/med",
               file.name="med.F")
