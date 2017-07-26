## Sample plots for low-dimensional simulation
library(bestsubset)
n = 100; p = 10
file.list = system(paste0("ls rds/sim.n",n,".p",p,".*.rds"),intern=TRUE)
method.nums = c(3,2,1,4)
method.names = c("Best subset","Forward stepwise","Lasso","Relaxed lasso")

beta.type = 2
rho = 0.35
short.list = grep(gsub("\\.","\\\\.",sprintf("*beta%i.rho%0.2f",beta.type,rho)),
                  file.list, val=TRUE)

plot.from.file(short.list, what="risk", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="none", make.pdf=TRUE, fig.dir="fig",
               file.name="lo.risk", h=4, w=4)

plot.from.file(short.list, what="error", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="none", make.pdf=TRUE, fig.dir="fig",
               file.name="lo.err", h=4, w=4)

plot.from.file(short.list, what="prop", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="none", make.pdf=TRUE, fig.dir="fig",
               file.name="lo.prop", h=4, w=4)

plot.from.file(short.list, what="nonzero", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="right", make.pdf=TRUE, fig.dir="fig",
               file.name="lo.nzs", h=4, w=5.5)
