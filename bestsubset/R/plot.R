#' Plot the results over several simulation settings.
#'
#' Plot the results over several sets of simulations, where the same methods are
#'   run over different simulations settings.
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
#' @param what one of "error" or "nonzeros", indicating whether to plot the
#'   relative test error or number of nonzeros, for each method across the given
#'   SNR levels. Default is "error".
#' @param tuning one of "validation" or "oracle", indicating whether the tuning
#'   parameter for each method should be chosen according to minimizing
#'   validation error, or according to minimizing test error. Default is
#'   "validation".
#' @param type Either "ave" or "med", indicating whether the average or median
#'   of the test error (or number of nonzeros, if what is "nonzeros") should be
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
#' @param pve Should the (population) proportion of variance explained be shown,
#'   corresponding to each SNR level under consideration? This is snr/(1+snr).
#'   Default is TRUE.
#' @param cols,main,cex.main,log,legend.pos graphical parameters.
#' @param tex.dir The latex directory to use. Default is ".".
#'
#' @export plot.many.sims

plot.many.sims = function(file.list, grouping, snr.vec,method.nums=NULL,
                          method.names=NULL, what=c("error","nonzeros"),
                          tuning=c("validation","oracle"), type=c("ave","med"),
                          std=TRUE, fig.dir=".", file.name=NULL, w=6, h=6,
                          mar=NULL, pve=TRUE, cols=1:8, main=NULL,
                          cex.main=1.25, log="x", legend.pos="bottomright",
                          tex.dir=".") {

  what = match.arg(what)
  tuning = match.arg(tuning)
  type = match.arg(type)
  if (is.null(file.name)) file.name = paste0("sim",1:length(file.list))
  if (is.null(mar)) {
    mar = c(4.25,4.25,1,1)
    if (pve) mar[3] = mar[3]+2
    if (!is.null(main)) mar[3] = mar[3]+2.25
  }
  if (is.null(main)) main = ""
  groups = unique(grouping)
  main = rep(main,length(groups))
  legend.pos = rep(legend.pos,length(groups))
  plot.name = numeric(length(groups))
  
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
      if (what=="error" && tuning=="validation" && type=="ave") {
        ymat = rbind(ymat, unlist(sim.obj$err.rel.tun.val.ave[method.nums]))
        ystd = rbind(ystd, unlist(sim.obj$err.rel.tun.val.std[method.nums]))
      }
      else if (what=="error" && tuning=="oracle" && type=="ave") {
        ymat = rbind(ymat, unlist(sim.obj$err.rel.tun.ora.ave[method.nums]))
        ystd = rbind(ystd, unlist(sim.obj$err.rel.tun.ora.std[method.nums]))
      }
      else if (what=="error" && tuning=="validation" && type=="med") {
        ymat = rbind(ymat, unlist(sim.obj$err.rel.tun.val.med[method.nums]))
        ystd = rbind(ystd, unlist(sim.obj$err.rel.tun.val.mad[method.nums]))
      }
      else if (what=="error" && tuning=="oracle" && type=="med") {
        ymat = rbind(ymat, unlist(sim.obj$err.rel.tun.ora.med[method.nums]))
        ystd = rbind(ystd, unlist(sim.obj$err.rel.tun.ora.mad[method.nums]))
      }
      else if (what=="nonzeros" && tuning=="validation" && type=="ave") {
        ymat = rbind(ymat, unlist(sim.obj$nzs.tun.val.ave[method.nums]))
        ystd = rbind(ystd, unlist(sim.obj$nzs.tun.val.std[method.nums]))
      }
      else if (what=="nonzeros" && tuning=="oracle" && type=="ave") {
        ymat = rbind(ymat, unlist(sim.obj$nzs.tun.ora.ave[method.nums]))
        ystd = rbind(ystd, unlist(sim.obj$nzs.tun.ora.std[method.nums]))
      }
      else if (what=="nonzeros" && tuning=="validation" && type=="med") {
        ymat = rbind(ymat, unlist(sim.obj$nzs.tun.val.med[method.nums]))
        ystd = rbind(ystd, unlist(sim.obj$nzs.tun.val.mad[method.nums]))
      }
      else if (what=="nonzeros" && tuning=="oracle" && type=="med") {
        ymat = rbind(ymat, unlist(sim.obj$nzs.tun.ora.med[method.nums]))
        ystd = rbind(ystd, unlist(sim.obj$nzs.tun.ora.mad[method.nums]))
      }
    }
      
    xlab = "Signal-to-noise ratio"
    ylab = ifelse(what=="error","(Test error)/sigma^2","Number of nonzeros")
    plot.name[i] = paste0(ifelse(what=="error","err.","nzs."),file.name[i])
    cols = rep(cols,length=length(method.nums))
    
    pdf(file=sprintf("%s/%s.pdf",fig.dir,plot.name[i]),w,h)
    par(mar=mar)
    matplot(snr.vec, ymat, type="o", pch=19, lty=1, col=cols, log=log,
            ylim=range(c(ymat-ystd,ymat+ystd)), xlab=xlab, ylab=ylab)
    for (k in 1:nrow(ystd)) {
      segments(snr.vec[k], ymat[k,]-ystd[k,], snr.vec[k],
               ymat[k,]+ystd[k,], lwd=0.5, col=cols)
    }
    
    if (legend.pos[i] != "") {
      legend(legend.pos[i],legend=method.names,col=cols,lty=1)
    }
    
    if (pve) {
      axis(side=3, at=snr.vec, labels=round(snr.vec/(1+snr.vec),2))
      mtext(side=3, line=2.5, "Proportion of variance explained")
    }

    mtext(side=3, line=4, main[i], cex=cex.main)
    graphics.off()
    cat(sprintf("%s/%s.pdf",fig.dir,plot.name[i]),"\n")
  }

  # Produce tex lines
  stem = get.stem(file.name)
  tex.name = paste0(ifelse(what=="error","err.","nzs."),
                    ifelse(tuning=="validation","val.","ora."),stem)
  tex.file = sprintf("%s/%s.tex",tex.dir,tex.name)
  cat("", file=tex.file, append=FALSE) # Wipe the file
  frac = 0.32
  for (i in 1:length(groups)) {
    cat(sprintf("\\includegraphics[width=%0.2f\\textwidth]{{%s/%s}.pdf}\n",
                frac,fig.dir,plot.name[i]), file=tex.file, append=TRUE)
  }
  cat(tex.file,"\n")
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
