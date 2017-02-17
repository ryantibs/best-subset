#' Plot the results over several simulation settings.
#'
#' Plot results over several sets of simulations, where the same methods are run
#'   over different simulations settings.
#'
#' @param file.list vector of strings that point to saved sim objects (each
#'   object produced by a call to \code{\link{sim.master}}).
#' @param grouping integer or factor vector indicating the grouping to use for
#'   the simulations. Within each group, the relative test error achieved by
#'   each method is plotted across the available SNR levels.
#' @param snr.vec Vector giving the SNR levels considered within each group.
#' @param method.nums the indices of the methods that should be plotted. Default
#'   is NULL, in which case all methods are plotted.
#' @param method.names the names of the methods that should be plotted. Default
#'   is NULL, in which case the names are extracted from the sim objects.
#' @param what one of "error", "prop", or "nonzero", indicating whether to plot
#'   the relative test error, test proportion of variance explained, or number
#'   of nonzeros, for each method across the given SNR levels. When what is
#'   "prop", the x-axis is population proportion of variance explained, instead
#'   of SNR. Default is "error". 
#' @param tuning one of "validation" or "oracle", indicating whether the tuning
#'   parameter for each method should be chosen according to minimizing
#'   validation error, or according to minimizing test error. Default is
#'   "validation".
#' @param type Either "ave" or "med", indicating whether the average or median
#'   of the test error (or number of nonzeros, if what is "nonzero") should be
#'   displayed. Default is "ave".
#' @param std Should standard errors be displayed (in parantheses)? When type
#'   is set to "med", the median absolute deviations are shown in place of the
#'   standard errors. Default is TRUE.
#' @param fig.dir The figure directory to use. Default is ".".
#' @param file.name Vector of strings, giving the file names to use for the 
#'   saved figures. Default is NULL, in which case the names "sim1", "sim2",
#'   etc. are used. (Extensions of "pdf" are always appended to the given file
#'   names.)
#' @param w,h the width and height (in inches) for the plots. Defaults are 6 for
#'   both.
#' @param mar the margins to use for the plots. Default is NULL, in which case
#'   the margins are set automatically (depending on whether not main is NULL).
#' @param cols,main,cex.main,log,legend.pos graphical parameters.
#' @param tex.dir The latex directory to use. Default is NULL, which means that
#'   no lines of tex (includegraphics statements, for the figures created) will
#'   be produced.
#'
#' @export plot.many.sims

plot.many.sims = function(file.list, grouping, snr.vec, method.nums=NULL,
                          method.names=NULL, what=c("error","prop","nonzero"),
                          tuning=c("validation","oracle"), type=c("ave","med"),
                          std=TRUE, fig.dir=".", file.name=NULL, w=6, h=6,
                          mar=NULL, cols=1:8, main=NULL, cex.main=1.25,
                          log=ifelse(what=="prop","","x"),
                          legend.pos="bottomright", tex.dir=NULL) {

  what = match.arg(what)
  tuning = match.arg(tuning)
  type = match.arg(type)
  if (is.null(file.name)) file.name = paste0("sim",1:length(file.list))
  if (is.null(mar)) {
    mar = c(4.25,4.25,1,1)
    #if (pve) mar[3] = mar[3]+2
    if (!is.null(main)) mar[3] = mar[3]+2.25
  }
  if (is.null(main)) main = ""
  groups = unique(grouping)
  main = rep(main,length(groups))
  legend.pos = rep(legend.pos,length(groups))
  plot.name = numeric(length(groups))

  # Set the x-axis, and axes labels
  if (what=="error") {
    xvar = snr.vec
    xlab = "Signal-to-noise ratio"
    ylab = "(Test error)/sigma^2"
  }
  else if (what=="prop") {
    xvar = snr.vec/(1+snr.vec)
    xlab = "Population proportion of var explained"
    ylab = "Test proportion of var explained"
  }
  else {
    xvar = snr.vec
    xlab = "Signal-to-noise ratio"
    ylab = "Number of nonzeros"
  }
  
  for (i in 1:length(groups)) {
    # Extract metrics from all the simulations in the current group
    jj = which(grouping==groups[i])
    ymat = c(); ystd = c()
    for (j in jj) {
      sim.obj = readRDS(file.list[j])
      if (is.null(method.nums)) method.nums = 1:length(sim.obj$err.train.ave)
      if (is.null(method.names)) method.names =
                                   names(sim.obj$err.train.ave[method.nums])

      # Grab the y-values for the plot
      comp.name = paste0(
        switch(what,error="err.rel.",prop="prop.",nonzero="nzs."),"tun.",
        switch(tuning,validation="val.",oracle="ora."))
      comp.name.val = paste0(comp.name, switch(type,ave="ave",med="med"))
      comp.name.std = paste0(comp.name, switch(type,ave="std",med="mad"))
      ymat = rbind(ymat, unlist(sim.obj[[comp.name.val]][method.nums]))
      ystd = rbind(ystd, unlist(sim.obj[[comp.name.std]][method.nums]))
    }

    # Set the plot name and the colors
    plot.name[i] = paste0(switch(what,error="err.",prop="pro.",nonzero="nzs."),
                          file.name[i])
    cols = rep(cols,length=length(method.nums))

    # Produce the plot
    pdf(file=sprintf("%s/%s.pdf",fig.dir,plot.name[i]),w,h)
    par(mar=mar)
    matplot(xvar, ymat, type="o", pch=19, lty=1, col=cols, log=log,
            ylim=range(c(ymat-ystd,ymat+ystd)), xlab=xlab, ylab=ylab,
            main=main[i], cex.main=cex.main)
    for (k in 1:nrow(ystd)) {
      segments(xvar[k], ymat[k,]-ystd[k,], xvar[k],
               ymat[k,]+ystd[k,], lwd=0.5, col=cols)
    }
    
    if (legend.pos[i] != "") {
      legend(legend.pos[i],legend=method.names,col=cols,lty=1)
    }
    
    ## if (pve) {
    ##   axis(side=3, at=xvar, labels=sprintf("%0.2f",snr.vec/(1+snr.vec)))
    ##   mtext(side=3, line=2.5, "Proportion of variance explained")
    ## }
    ## mtext(side=3, line=4, main[i], cex=cex.main)
    graphics.off()
    cat(sprintf("%s/%s.pdf",fig.dir,plot.name[i]),"\n")
  }

  # Produce tex lines, if we are asked to
  if (!is.null(tex.dir)) {
    stem = get.stem(file.name)
    tex.name = paste0(switch(what,error="err.",prop="pro.",nonzero="nzs."),
                      switch(tuning,validation="val.",,oracle="ora."),stem)
    tex.file = sprintf("%s/%s.tex",tex.dir,tex.name)
    cat("", file=tex.file, append=FALSE) # Wipe the file
    frac = 0.32
    for (i in 1:length(groups)) {
      cat(sprintf("\\includegraphics[width=%0.2f\\textwidth]{{%s/%s}.pdf}\n",
                  frac,fig.dir,plot.name[i]), file=tex.file, append=TRUE)
    }
    cat(tex.file,"\n")
  }
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
