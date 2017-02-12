#' Predictors and responses generation.
#'
#' Generate a predictor matrix x, and response vector y, following a specified
#'   setup.  Actually, three pairs of predictors and responses are generated:
#'   one for training, one for validation, and one for testing.
#'
#' @param n,p The number of training observations, and the number of predictors.
#' @param nval,ntest The number of validation observations, and the number of
#'   testing observations.
#' @param rho Parameter that drives pairwise correlations of the predictor
#'   variables; specifically, predictors i and j have population correlation
#'   rho^abs(i-j). Default is 0.
#' @param s number of nonzero coefficients in the underlying regression model.
#'   Default is 5. (Ignored if beta.type is 4, in which case the number of
#'   nonzero coefficients is 6; and if beta.type is 5, it is interpreted as a
#'   the number of strongly nonzero coefficients in a weak sparsity model.)
#' @param beta.type Integer taking values in between 1 and 5, used to specify
#'   the pattern of nonzero coefficients in the underlying regression model; see
#'   details below. Default is 1.
#' @param snr Desired signal-to-noise ratio (SNR), i.e., var(mu)/sigma^2 where
#'   mu is mean and sigma^2 is the error variance. The error variance is set so
#'   that the given SNR is achieved. Default is 1.
#'
#' @return A list with the following components: x, y, xval, yval, xtest, ytest,
#'   mutest, beta, and sigma.
#' 
#' @details The predictors are normal with covariance sigma^2 * Sigma, where
#'   sigma^2 is set according to the desired signal-to-noise ratio, and Sigma
#'   has (i,j)th entry rho^abs(i-j). The first 4 options for the nonzero pattern
#'   of the underlying regression coefficients beta follow the simulation setup
#'   in Bertsimas, King, and Mazumder (2016), and the last is a weak sparsity
#'   option:
#'   \itemize{
#'   \item 1: beta has s components of 1, occurring at (roughly) equally-spaced
#'      indices in between 1 and p
#'   \item 2: beta has its first s components equal to 1
#'   \item 3: beta has its first s components taking nonzero values, where the
#'       decay in a linear fashion from 10 to 0.5
#'   \item 4: beta has its first 6 components taking the nonzero values -10,-6,
#'       -2,2,6,10
#'   \item 5: beta has its first s components equal to 1, and the rest decaying
#'       to zero at an exponential rate
#'   }
#'
#' @author Trevor Hastie, Rob Tibshirani, Ryan Tibshirani
#' @references Simulation setup based on "Best subset selection via a modern
#'   optimization lens" by Dimitris Bertsimas, Angela King, and Rahul Mazumder,
#'   Annals of Statistics, 44(2), 813-852, 2016.
#' @example examples/ex.fs.R
#' @export sim.xy

sim.xy = function(n, p, nval, ntest, rho=0, s=5, beta.type=1, snr=1) {
  # Generate predictors
  x = matrix(rnorm(n*p),n,p)
  xval = matrix(rnorm(nval*p),nval,p)
  xtest = matrix(rnorm(ntest*p),ntest,p)

  # Introduce autocorrelation, if needed
  if (rho != 0) {
    inds = 1:p
    Sigma = rho^abs(outer(inds, inds, "-"))
    obj = svd(Sigma)
    Sigma.half = obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v)
    x = x %*% Sigma.half
    xval = xval %*% Sigma.half
    xtest = xtest %*% Sigma.half
  }

  # Generate underlying coefficients
  s = min(s,p)
  beta = rep(0,p)
  if (beta.type==1) {
    beta[round(seq(1,p,length=s))] = 1
  } else if (beta.type==2) {
    beta[1:s] = 1
  } else if (beta.type==3) {
    beta[1:s] = seq(10,0.5,length=s)
  } else if (beta.type==4) {
    beta[1:6] = c(-10,-6,-2,2,6,10)
  } else {
    beta[1:s] = 1
    beta[(s+1):p] = 0.5^(1:(p-s))
  }

  # Set snr based on sample variance on large test set 
  mutest = as.numeric(xtest %*% beta)
  vmu = var(mutest) 
  sigma = sqrt(vmu/snr)

  # Generate responses
  y = as.numeric(x %*% beta + rnorm(n)*sigma)
  yval = as.numeric(xval %*% beta + rnorm(nval)*sigma)
  ytest = as.numeric(mutest + rnorm(ntest)*sigma)
 
  enlist(x,y,xval,yval,xtest,ytest,mutest,beta,sigma)
}

#' Master function for running simulations.
#'
#' Run a set of simulations with the specified configuration.
#'
#' @param n,p The number of training observations, and the number of predictors.
#' @param nval,ntest The number of validation observations, and the number of
#'   testing observations.
#' @param reg.funs This is a list of functions, representing the regression
#'   procedures to be used (evaluated) in the simulation. Each element of the 
#'   list must be a function that takes x, y (the training predictor matrix and
#'   response vector) as its only two (mandatory) arguments, and must return an 
#'   object with associated coef and predict methods. The coef method must take
#'   obj (the returned object) and return a matrix of coefficients, with one
#'   column per tuning parameter value inherent to the regression method. The
#'   predict method must take obj, newx (the returned object and a new predictor
#'   matrix) and return a matrix of predictions, again with one column per
#'   tuning parameter value inherent to the regression method.
#' @param Number of repetitions of which to average the results. Default is 50.
#' @param seed Seed to be set for the overall random number generation, i.e.,
#'   set before repetitions are begun (for reproducibility of the simulation
#'   results). Default is NULL, which effectively sets no seed.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#' @param file,file.rep Name of a file to which simulation results are saved
#'   (using saveRDS), and a number of repetitions after which intermediate
#'   results are saved. Setting file to NULL is interpreted to mean that no
#'   simulations results should be saved; setting file.rep to 0 is interpreted
#'   to mean that simulations results should be saved at the very end, i.e., no
#'   intermediate saving. Defaults are NULL and 5, respectively.
#' @param rho,s,beta.type,snr Arguments to pass to \code{\link{sim.xy}}; see the
#'   latter's help file for details.
#'
#' @return A list with components err.train, err.val, err.test, err.rel, risk,
#'   nzs, opt for the training error, validation error, test error, relative
#'   test error (test error divided by sigma^2), risk, number of
#'   selected nonzero coefficients, and optimism (difference in test and
#'   training errors). These are each lists of length N, where N is the number
#'   of regression methods under consideration (the length of reg.funs). The
#'   ith element of each list is then a matrix of dimension nrep x m, where m 
#'   the number of tuning parameters inherent to the ith method. The returned
#'   components err.train.ave, err.val.ave, err.test.ave, err.rel.ave, risk.ave,
#'   nzs.ave, opt.ave return the averages of the training error, validation
#'   error, etc. over the nrep repetitions. Similarly for the components with
#'   postfixes .std, .med, .mad, which return the standard deviation, median,
#'   and median absolute deviation, respectively. The returned components with
#'   postfixes .tun.val and .tun.orc are matrices of dimension N x nrep, which
#'   return the training error, validation error, etc. when the tuning parameter
#'   for each regression method in each repetition is chosen by validation
#'   tuning (best validation error) or by oracle tuning (best test error).
#'
#' @seealso \code{\link{sim.xy}}
#' @author Trevor Hastie, Robert Tibshirani, Ryan Tibshirani
#' @references The structure of this simulation code based on that from the
#'   \code{conformalInference} package.
#' @example examples/ex.sim.master.R
#' @export sim.master

sim.master = function(n, p, nval, ntest, reg.funs, nrep=50, seed=NULL,
  verbose=FALSE, file=NULL, file.rep=5, rho=0, s=5, beta.type=1, snr=1) {
  
  this.call = match.call()
  if (!is.null(seed)) set.seed(seed)

  N = length(reg.funs)
  reg.names = names(reg.funs)
  if (is.null(reg.names)) reg.names = paste("Method",1:N)
  
  err.train = err.val = err.test = err.rel = risk = nzs = opt = runtime =
    vector(mode="list",length=N)
  names(err.train) = names(err.val) = names(err.test) = names(err.rel) = 
    names(risk) = names(nzs) = names(opt) = names(runtime) = reg.names
  for (j in 1:N) {
    err.train[[j]] = err.val[[j]] = err.test[[j]] = err.rel[[j]] =
      risk[[j]] = nzs[[j]] = opt[[j]] = runtime[[j]] = matrix(NA,nrep,1)
  }
  filled = rep(FALSE,N)

  # Loop through the repetitions
  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Simulation %i (of %i) ...\n",i,nrep))
      cat("  Generating data ...\n")
    }
    
    # Generate x, y, xval, yval, xtest, ytest
    xy.obj = sim.xy(n,p,nval,ntest,rho,s,beta.type,snr)

    # Loop through the regression methods
    for (j in 1:N) {
      if (verbose) {
        cat(sprintf("  Applying regression method %i (of %i) ...\n",
                    j,N))
      }
      
      tryCatch({
        # Apply the regression method in hand
        runtime[[j]][i] = system.time({   
          reg.obj = reg.funs[[j]](xy.obj$x,xy.obj$y)
        })[1]

        # Grab the coefficients, and predicted values on the training,
        # validation, and testing sets
        beta = as.matrix(coef(reg.obj))
        muhat.train = as.matrix(predict(reg.obj,xy.obj$x))
        muhat.val = as.matrix(predict(reg.obj,xy.obj$xval))
        muhat.test = as.matrix(predict(reg.obj,xy.obj$xtest))
        
        # Populate empty matrices for our metrics, of appropriate dimension
        if (!filled[j]) {
          err.train[[j]] = err.val[[j]] = err.test[[j]] = err.rel[[j]] =
            risk[[j]] = nzs[[j]] = opt[[j]] = matrix(NA,nrep,ncol(beta))
          filled[j] = TRUE
          # N.B. Filling with NAs is important, because the filled flag could
          # be false for two reasons: i) we are at the first iteration, or ii)
          # we've failed in all previous iters to run the regression method
        }

        # Record all of our metrics
        err.train[[j]][i,] = colMeans((muhat.train - xy.obj$y)^2)
        err.val[[j]][i,] = colMeans((muhat.val - xy.obj$yval)^2)
        err.test[[j]][i,] = colMeans((muhat.test - xy.obj$ytest)^2)
        err.rel[[j]][i,] = err.test[[j]][i,] / xy.obj$sigma^2
        risk[[j]][i,] = colMeans((muhat.test - xy.obj$mutest)^2)
        nzs[[j]][i,] = colSums(beta != 0)
        opt[[j]][i,] = err.test[[j]][i,] - err.train[[j]][i,]
      }, error = function(err) {
        if (verbose) {
          cat(paste("    Oops! Something went wrong, see error message",
                    "below; recording all metrics here as NAs ...\n"))
          cat("    ***** Error message *****\n")
          cat(sprintf("    %s\n",err$message))
          cat("    *** End error message ***\n")
        }
        # N.B. No need to do anything, the metrics are already filled with NAs
      })
    }

    # Save intermediate results?
    if (!is.null(file) && file.rep > 0 && i %% file.rep == 0) {
      saveRDS(enlist(err.train,err.val,err.test,err.rel,risk,nzs,
                     opt,runtime),file)
    }
  }

  # Save results now (in case of an error that might occur below)
  out = enlist(err.train,err.val,err.test,err.rel,risk,nzs,opt,runtime)
  if (!is.null(file)) saveRDS(out, file)
  
  # Tune according to validation error, and according to test error
  out = tune.methods(out)

  # Aggregate our metrics over the simulations
  out = compute.stats(out)
  
  # Save final results
  out = c(out,list(rho=rho,s=s,beta.type=beta.type,snr=snr,call=this.call))
  class(out) = "sim"
  if (!is.null(file)) { saveRDS(out, file); invisible(out) }
  else return(out)
}

##############################

tune.methods = function(obj) {
  # Tune first by min validation error
  N = length(obj$err.val) # Number of methods
  nrep = nrow(obj$err.val[[1]]) # Number of repetitions
  method.names = names(obj$err.val) # Names of methods
  tun.val = err.rel.tun.val = nzs.tun.val = vector(mode="list",length=N)
  names(tun.val) = names(err.rel.tun.val) = names(nzs.tun.val) = method.names
  for (j in 1:N) {
    tun.val[[j]] = err.rel.tun.val[[j]] = nzs.tun.val[[j]] = matrix(NA,nrep,1)
  }
  
  for (i in 1:nrep) {
    for (j in 1:N) {
      tun.val[[j]][i] = which.min(obj$err.val[[j]][i,])
      if (length(tun.val[[j]][i]) == 0) { # Error occurred in rep i
        tun.val[[j]][i] = NA
        err.rel.tun.val[[j]][i] = NA
        nzs.tun.val[[j]][i] = NA
      }
      else { # Repetition i was completed as usual
        err.rel.tun.val[[j]][i] = obj$err.rel[[j]][i,tun.val[[j]][i]]
        nzs.tun.val[[j]][i] = obj$nzs[[j]][i,tun.val[[j]][i]]
      }
    }
  }

  # Tune second by min test error
  tun.ora = err.rel.tun.ora = nzs.tun.ora = vector(mode="list",length=N)
  names(tun.ora) = names(err.rel.tun.ora) = names(nzs.tun.ora) = method.names
  for (j in 1:N) {
    tun.ora[[j]] = err.rel.tun.ora[[j]] = nzs.tun.ora[[j]] = matrix(NA,nrep,1)
  }
  
  for (j in 1:N) {
    for (i in 1:nrep) {
      tun.ora[[j]][i] = which.min(obj$err.test[[j]][i,])
      if (length(tun.ora[[j]][i]) == 0) { # Error occurred in rep i
        tun.ora[[j]][i] = NA
        err.rel.tun.ora[[j]][i] = NA
        nzs.tun.ora[[j]][i] = NA
      }
      else { # Repetition i was completed as usual
        err.rel.tun.ora[[j]][i] = obj$err.rel[[j]][i,tun.ora[[j]][i]]
        nzs.tun.ora[[j]][i] = obj$nzs[[j]][i,tun.ora[[j]][i]]
      }
    }
  }

  return(c(obj,enlist(tun.val,err.rel.tun.val,nzs.tun.val,
                      tun.ora,err.rel.tun.ora,nzs.tun.ora)))
}

compute.stats = function(obj) {
  nmet = length(obj) # Number of metrics
  N = length(obj$err.val) # Number of methods
  obj.ave = obj.std = vector(mode="list", length=nmet)
  obj.med = obj.mad = vector(mode="list", length=nmet)
  names(obj.ave) = paste0(names(obj), ".ave")
  names(obj.std) = paste0(names(obj), ".std")
  names(obj.med) = paste0(names(obj), ".med")
  names(obj.mad) = paste0(names(obj), ".mad")
  
  for (i in 1:nmet) {
    obj.ave[[i]] = obj.std[[i]] = obj.med[[i]] = obj.mad[[i]] =
      vector(mode="list", length=N)
    names(obj.ave[[i]]) = names(obj.std[[i]]) = names(obj.med[[i]]) =
      names(obj.mad[[i]]) = names(obj[[i]])
    
    for (j in 1:N) {
      obj.ave[[i]][[j]] = colMeans(obj[[i]][[j]], na.rm=TRUE)
      obj.std[[i]][[j]] = apply(obj[[i]][[j]], 2, sd, na.rm=TRUE) /
        sqrt(colSums(!is.na(obj[[i]][[j]])))
      obj.med[[i]][[j]] = apply(obj[[i]][[j]], 2, median, na.rm=TRUE)
      obj.mad[[i]][[j]] = apply(obj[[i]][[j]], 2, mad, na.rm=TRUE) /
        sqrt(colSums(!is.na(obj[[i]][[j]])))
    }
  }
  
  return(c(obj,obj.ave,obj.std,obj.med,obj.mad))
}

##############################

#' Print function for sim object.
#'
#' Summarize and print the results of a set of simulations, stored an object
#'   of class sim (produced by \code{\link{sim.master}}).
#'
#' @param x The sim object.
#' @param type Either "ave" or "med", indicating whether the average or median
#'   of the relative test error metric should be displayed. Default is "ave".
#' @param std Should standard errors be displayed (in parantheses)? When type
#'   is set to "med", the median absolute deviations are shown in place of the
#'   standard errors. Default is TRUE. 
#' @param digits Number of digits to display. Default is 3. 
#' @param ... Other arguments (currently not used).
#'
#' @export print.sim
#' @export 

print.sim = function(x, type=c("ave","med"), std=TRUE, digits=3, ...) {
  type = match.arg(type)
  
  if (!is.null(x$call)) {
    cat("\nCall:\n") # TODO append call
    dput(x$call)
  }
  N = length(x$err.rel.tun.val.ave)
  nrep = length(x$err.rel.tun.val.ave[[1]])

  cat("\nResults for tuning parameters chosen based on validation set:\n\n")
  if (type=="ave") {
    col1 = unlist(x$err.rel.tun.val.ave)
    col2 = unlist(x$nzs.tun.val.ave)
    col1.std = unlist(x$err.rel.tun.val.std)
    col2.std = unlist(x$nzs.tun.val.std)
  }
  else {
    col1 = unlist(x$err.rel.tun.val.med)
    col2 = unlist(x$nzs.tun.val.med)
    col1.std = unlist(x$err.rel.tun.val.mad)
    col2.std = unlist(x$nzs.tun.val.mad)
  }

  tab = round(cbind(col1,col2),digits)
  tab.std = round(cbind(col1.std,col2.std), digits)
  if (std) tab = matrix(paste0(tab," (",tab.std,")"),ncol=2)
  rownames(tab) = names(x$err.rel.tun.val.ave)
  colnames(tab) = c("(Test error)/sigma^2","Nonzero coefficients")
  print(tab,quote=F)

  cat("\nResults for tuning parameters chosen based on test set (oracle):\n\n")
  if (type=="ave") {
    col1 = unlist(x$err.rel.tun.ora.ave)
    col2 = unlist(x$nzs.tun.ora.ave)
    col1.std = unlist(x$err.rel.tun.ora.std)
    col2.std = unlist(x$nzs.tun.ora.std)
  }
  else {
    col1 = unlist(x$err.rel.tun.ora.med)
    col2 = unlist(x$nzs.tun.ora.med)
    col1.std = unlist(x$err.rel.tun.ora.mad)
    col2.std = unlist(x$nzs.tun.ora.mad)
  }

  tab = round(cbind(col1,col2),digits)
  tab.std = round(cbind(col1.std,col2.std), digits)
  if (std) tab = matrix(paste0(tab," (",tab.std,")"),ncol=2)
  rownames(tab) = names(x$err.rel.tun.ora.ave)
  colnames(tab) = c("(Test error)/sigma^2","Nonzero coefficients")
  print(tab,quote=F)

  cat("\n")
  invisible()
}

#' Plot function for sim object.
#'
#' Plot the results of a set of simulations, stored in an object of class sim 
#'   (produced by \code{\link{sim.master}}).
#'
#' @param x The sim object.
#' @param method.nums the indices of the methods that should be plotted. Default
#'   is to 1:length(x$err.rel.ave), which plots all methods.
#' @param method.names the names of the methods that should be plotted. Default
#'   is NULL, in which case the names are extracted from the sim object.
#' @param type Either "ave" or "med", indicating whether the average or median
#'   of the relative test error metric should be displayed. Default is "ave".
#' @param std Should standard errors be displayed (in parantheses)? When type
#'   is set to "med", the median absolute deviations are shown in place of the
#'   standard errors. Default is TRUE.
#' @param cols,main,cex.main,legend.pos graphical parameters.
#' @param make.pdf Should a pdf be produced? Default is FALSE.
#' @param fig.dir,file.name The figure directory and file name to use, only 
#'   when make.pdf is TRUE. Defaults are "." and "sim". (An extension of "pdf"
#'   is always appended to the given file name.)
#' @param w,h the width and height (in inches) for the plot, used only when
#'   make.pdf is TRUE. Defaults are 6 for both.
#' @param mar the margins to use for the plot. Default is NULL, in which case
#'   the margins are set automatically (depending on whether not main is NULL).
#'
#' @export plot.sim
#' @export

plot.sim = function(x, method.nums=1:length(x$err.train.ave), method.names=NULL,
                    type=c("ave","med"), std=TRUE, cols=1:8, main=NULL,
                    cex.main=1.25, legend.pos=c("topright"),
                    make.pdf=FALSE, fig.dir=".", file.name="sim", w=6, h=6, 
                    mar=NULL) {

  type = match.arg(type)
  if (is.null(method.names)) method.names = names(x$err.train.ave[method.nums])
  if (is.null(mar) && is.null(main)) mar = c(4.25,4.25,1,1)
  if (is.null(mar) && !is.null(main)) mar = c(4.25,4.25,2.25,1)
  if (is.null(main)) main = ""
  cols = rep(cols,length=length(method.nums))
  ii = method.nums

  if (type=="ave") {
    xlist = x$nzs.ave
    ylist = x$err.rel.ave
    ybars = x$err.rel.std
  }
  else {
    xlist = x$nzs.med
    ylist = x$err.rel.med
    ybars = x$err.rel.mad
  }
  
  xlab = "Number of nonzero coefficients"
  ylab = "(Test error)/sigma^2"
  
  if (make.pdf) pdf(file=sprintf("%s/%s.pdf",fig.dir,file.name),w,h)
  par(mar=mar)
  xlim = range(unlist(xlist),na.rm=TRUE)
  ylim = range(unlist(ylist),na.rm=TRUE)
  if (std) ylim = range(c(unlist(ylist)-unlist(ybars),
                          unlist(ylist)+unlist(ybars)),na.rm=TRUE)
  plot(xlist[[ii[1]]], ylist[[ii[1]]], col=cols[1], type="o",
       xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main,
       cex.main=cex.main)
  if (std) segments(xlist[[ii[1]]], ylist[[ii[1]]]-ybars[[ii[1]]],
                    xlist[[ii[1]]], ylist[[ii[1]]]+ybars[[ii[1]]],
                    lwd=0.5, col=cols[1])   
  for (i in Seq(2,length(ii))) {
    points(xlist[[ii[i]]], ylist[[ii[i]]], col=cols[i], type="o")
    if (std) segments(xlist[[ii[i]]], ylist[[ii[i]]]-ybars[[ii[i]]],
                      xlist[[ii[i]]], ylist[[ii[i]]]+ybars[[ii[i]]],
                      lwd=0.5, col=cols[i])
  }
  if (legend.pos != "") legend(legend.pos,legend=method.names,col=cols,lty=1)
  if (make.pdf) graphics.off()
}

#' Print function for latex-style tables.
#'
#' Print a given table in format digestable by latex.
#'
#' @export print.tex

print.tex = function(tab, tab.se=NULL, digits=3, file=NULL, align="l") {
  tab = round(tab,digits)
  n = nrow(tab); m = ncol(tab)
  rownms = rownames(tab)
  colnms = colnames(tab)
  if (!is.null(tab.se)) {
    tab = matrix(paste0(tab, " (", round(tab.se,digits), ")"), ncol=m)
  }
  
  if (is.null(file)) file = ""
  cat("", file=file, append=F) # Wipe the file, or the console
  cat(paste0("\\begin{tabular}{|",paste0(rep(align,m+1),collapse="|"),"|}\n"),
      file=file, append=T)
  cat("\\hline\n", file=file, append=T)
  if (!is.null(colnms)) {
    for (j in 1:m) cat(paste0("& ", colnms[j]," "), file=file, append=T)
    cat("\\\\\n", file=file, append=T)
    cat("\\hline\n", file=file, append=T)
  }
  for (i in 1:n) {
    cat(paste0(rownms[i], " "), file=file, append=T)
    for (j in 1:m) cat(paste0("& ",tab[i,j]," "), file=file, append=T)
    cat("\\\\\n", file=file, append=T)
    cat("\\hline\n", file=file, append=T)
  }
  cat(paste0("\\end{tabular}\n"), file=file, append=T)
}
