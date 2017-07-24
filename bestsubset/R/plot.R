#' Plot the results over several simulation settings.
#'
#' Plot results over several sets of simulations saved in files, where the same
#'   methods are run over different simulations settings.
#'
#' @param file.list Vector of strings that point to saved sim objects (each
#'   object produced by a call to \code{\link{sim.master}}).
#' @param row,col One of "beta", "rho", or "snr", indicating the variables for
#'   the rows and columns of plotting grid; note that row and col must be
#'   different.  Default is row="beta" and col="rho", so that the plotting grid
#'   displays the metric specified in what (default is relative test error, see
#'   below) versus the SNR, in a plotting grid with the coefficient types across
#'   the rows, and the correlation levels across the columns.
#' @param method.nums The indices of the methods that should be plotted. Default
#'   is NULL, in which case all methods are plotted.
#' @param method.names The names of the methods that should be plotted. Default
#'   is NULL, in which case the names are extracted from the sim objects.
#' @param what One of "error", "risk", "prop", or "nonzero", indicating whether to plot
#'   the relative test error, test proportion of variance explained, or number
#'   of nonzeros, for each method across the given SNR levels. When what is
#'   "prop", the x-axis is population proportion of variance explained, instead
#'   of SNR. Default is "error". 
#' @param tuning Either "validation" or "oracle", indicating whether the tuning
#'   parameter for each method should be chosen according to minimizing
#'   validation error, or according to minimizing test error. Default is
#'   "validation".
#' @param type Either "ave" or "med", indicating whether the average or median
#'   of the test error (or number of nonzeros, if what is "nonzero") should be
#'   displayed. Default is "ave".
#' @param std Should standard errors be displayed (in parantheses)? When type
#'   is set to "med", the median absolute deviations are shown in place of the
#'   standard errors. Default is TRUE.
#' @param lwd,pch,main graphical parameters.
#' @param make.pdf Should a pdf be produced? Default is FALSE.
#' @param fig.dir,file.name The figure directory and file name to use, only
#'   when make.pdf is TRUE. Defaults are "." and "sim". (An extension of "pdf"
#'   is always appended to the given file name.)
#' @param w,h the width and height (in inches) for the plot, used only when
#'   make.pdf is TRUE. Defaults are 6 for both.
#'
#' @export plot.from.file

plot.from.file = function(file.list,
                          row=c("beta","rho","snr"), col=c("rho","beta","snr"),
                          method.nums=NULL, method.names=NULL,
                          what=c("error","risk","prop","nonzero"), rel.to=NULL,
                          tuning=c("validation","oracle"), type=c("ave","med"),
                          std=TRUE, lwd=1, pch=19, main=NULL, make.pdf=FALSE,
                          fig.dir=".", file.name="sim", w=8, h=10) {

  # Check for ggplot2 package
  if (!require("ggplot2",quietly=TRUE)) {
    stop("Package ggplot2 not installed (required here)!")
  }
  
  row = match.arg(row)
  col = match.arg(col)
  if (row==col) stop("row and col must be different")
  
  what = match.arg(what)
  tuning = match.arg(tuning)
  type = match.arg(type)

  # Set the method numbers and names
  sim.obj = readRDS(file.list[1])
  if (is.null(method.nums)) method.nums = 1:length(sim.obj$err.test)
  if (is.null(method.names)) method.names =
                               names(sim.obj$err.test[method.nums])
  N = length(method.nums)
  
  # Set the base number and name
  if (is.null(rel.to)) {
    base.num = 0
    base.name = "null"
  }
  else {
    base.num = which(method.nums==rel.to)
    base.name = tolower(method.names[base.num])
  }
  
  # Set the y-label
  ylab = switch(what,
                error=paste0("Relative test error (to ",base.name,")"),
                risk=paste0("Relative risk (to ",base.name,")"),
                prop="Proportion of variance explained",
                nonzero="Number of nonzeros")

  # Collect the y-variable from the file list
  yvec = ybar = beta.vec = rho.vec = snr.vec = c()
  for (i in 1:length(file.list)) {
    sim.obj = readRDS(file.list[i])
    beta.vec = c(beta.vec,rep(sim.obj$beta.type,N))
    rho.vec = c(rho.vec,rep(sim.obj$rho,N))
    snr.vec = c(snr.vec,rep(sim.obj$snr,N))

    z = sim.obj[[switch(what,
                        error="err.test",
                        risk="risk",
                        prop="prop",
                        nonzero="nzs")]]
    res = tune.and.aggregate(sim.obj, z)

    # For prop and nonzero we ignore any request for a relative metric
    if (what=="prop" || what=="nonzero") {
      yvec = c(yvec,res[[paste0("z.",substr(tuning,1,3),".",type)]][method.nums])
      ybar = c(ybar,res[[paste0("z.",substr(tuning,1,3),".",
                                ifelse(type=="ave","std","mad"))]][method.nums])
    }

    # For err and risk we respect the request for a relative metric
    else {
      # First build the relative metric
      met = res[[paste0("z.",substr(tuning,1,3))]][method.nums]
      if (base.num == 0 && what=="error") denom = sim.obj$err.null
      else if (base.num == 0 && what=="risk") denom = sim.obj$risk.null
      else denom = met[[base.num]]
      z.rel = lapply(met, function(v) v / denom)
      # Now aggregate the relative metric
      res2 = tune.and.aggregate(sim.obj, z.rel, tune=FALSE)
      yvec = c(yvec,unlist(res2[[paste0("z.",type)]]))
      ybar = c(ybar,unlist(res2[[paste0("z.",ifelse(type=="ave",
                                                        "std","mad"))]]))
    }
  }
  
  # Set the x-variable and x-label
  xvec = snr.vec
  xlab = "Signal-to-noise ratio"

  # Produce the plot
  beta.vec = factor(beta.vec)
  rho.vec = factor(rho.vec)
  snr.vec = factor(snr.vec)
  levels(beta.vec) = paste("Beta-type", levels(beta.vec))
  levels(rho.vec) = paste("Correlation", levels(rho.vec))
  
  dat = data.frame(x=xvec, y=yvec, se=ybar,
                   beta=beta.vec, rho=rho.vec, snr=snr.vec,
                   Method=factor(rep(method.names, length=length(xvec))))

  gp = ggplot(dat, aes(x=xvec,y=yvec,color=Method)) +
    xlab(xlab) + ylab(ylab) +
    geom_line(lwd=lwd) + geom_point(pch=pch) +
    facet_grid(formula(paste(row,"~",col))) +
    theme_bw() + theme(legend.pos="bottom")
  if (!("snr" %in% c(row,col))) {
    # If SNR is being plotted on the x-axis in each plot, then define special
    # x-axis ticks and put the x-axis on a log scale
    snr.breaks = round(exp(seq(from=min(log(xvec)),
                               to=max(log(xvec)),length=4)),2)
    gp = gp + scale_x_continuous(trans="log", breaks=snr.breaks)
  }
  if (std) gp = gp + geom_errorbar(aes(ymin=yvec-se,ymax=yvec+se), width=0.02)
  if (what=="prop") gp = gp + geom_line(aes(x=xvec, y=xvec/(1+xvec)), lwd=0.5,
                                        linetype=3, color="black")
  if (what =="nonzero") gp = gp + geom_line(aes(x=xvec, y=sim.obj$s), lwd=0.5,
                                            linetype=3, color="black")
  if (!is.null(main)) gp = gp + ggtitle(main) 
  if (make.pdf) ggsave(sprintf("%s/%s.pdf",fig.dir,file.name),
                       height=h, width=w, device="pdf")
  else gp
}

get.stem = function(str.vec) {
  str.list = strsplit(str.vec, split="\\.")
  k = 0
  while (TRUE) {
    vec = c()
    for (i in 1:length(str.list)) {
      if (length(str.list[[i]]) < k+1) break
      vec = c(vec, str.list[[i]][k+1])
    }
    if (length(vec) < length(str.list)) break
    if (length(unique(vec)) > 1) break
    k = k+1
  }
  if (k == 0) stem = "foo"
  else stem = paste0(str.list[[1]][1:k], collapse=".")
  return(stem)
}
