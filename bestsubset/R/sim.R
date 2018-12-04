#' Predictors and responses generation.
#'
#' Generate a predictor matrix x, and response vector y, following a specified
#'   setup.  Actually, two pairs of predictors and responses are generated:
#'   one for training, and one for validation.
#'
#' @param n,p The number of training observations, and the number of predictors.
#' @param nval The number of validation observations.
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
#' @return A list with the following components: x, y, xval, yval, Sigma, beta,
#'   and sigma.
#'
#' @details The data model is: \eqn{Y \sim N(X\beta, \sigma^2 I)}.
#'   The predictor variables have covariance matrix Sigma, with (i,j)th entry
#'   rho^abs(i-j). The error variance sigma^2 is set according to the desired
#'   signal-to-noise ratio. The first 4 options for the nonzero pattern
#'   of the underlying regression coefficients beta follow the simulation setup
#'   in Bertsimas, King, and Mazumder (2016), and the 5th is a weak sparsity
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

sim.xy = function(n, p, nval, rho=0, s=5, beta.type=1, snr=1) {
  # Generate predictors
  x = matrix(rnorm(n*p),n,p)
  xval = matrix(rnorm(nval*p),nval,p)

  # Introduce autocorrelation, if needed
  if (rho != 0) {
    inds = 1:p
    Sigma = rho^abs(outer(inds, inds, "-"))
    obj = svd(Sigma)
    Sigma.half = obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v)
    x = x %*% Sigma.half
    xval = xval %*% Sigma.half
  }
  else Sigma = diag(1,p)

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

  # Set snr based on sample variance on infinitely large test set
  vmu = as.numeric(t(beta) %*% Sigma %*% beta)
  sigma = sqrt(vmu/snr)

  # Generate responses
  y = as.numeric(x %*% beta + rnorm(n)*sigma)
  yval = as.numeric(xval %*% beta + rnorm(nval)*sigma)

  enlist(x,y,xval,yval,Sigma,beta,sigma)
}

#' Master function for running simulations.
#'
#' Run a set of simulations with the specified configuration.
#'
#' @param n,p The number of training observations, and the number of predictors.
#' @param nval The number of validation observations.
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
#' @param rho,s,beta.type,snr. Arguments to pass to \code{\link{sim.xy}}; see
#'   the latter's help file for details.
#'
#' @return A list with components err.train, err.val, err.test, prop, risk, nzs,
#'   fpos, fneg, F1, opt for the training error, validation error, test error,
#'   test proportion of variance explained, risk, number of selected nonzero
#'   coefficients, number of false positives, number of false negatives, F1
#'   measure, and relative optimism (difference in test error and training
#'   error, divided by training error), respectively.  These are each lists of
#'   length N, where N is the number of regression methods under consideration
#'   (the length of reg.funs). The ith element of each list is then a matrix of
#'   dimension nrep x m, where m is the number of tuning parameters inherent to
#'   the ith method.
#'
#' @seealso \code{\link{sim.xy}}
#' @author Trevor Hastie, Robert Tibshirani, Ryan Tibshirani
#' @references The structure of this simulation code based on that from the
#'   \code{conformalInference} package.
#' @example examples/ex.sim.master.R
#' @export sim.master

sim.master = function(n, p, nval, reg.funs, nrep=50, seed=NULL, verbose=FALSE,
                      file=NULL, file.rep=5, rho=0, s=5, beta.type=1, snr=1) {
  
  this.call = match.call()
  if (!is.null(seed)) set.seed(seed)

  N = length(reg.funs)
  reg.names = names(reg.funs)
  if (is.null(reg.names)) reg.names = paste("Method",1:N)

  err.train = err.val = err.test = prop = risk = nzs = fpos = fneg = F1 = 
    opt = runtime = vector(mode="list",length=N)
  names(err.train) = names(err.val) = names(err.test) = names(prop) =
    names(risk) = names(nzs) = names(fpos) = names(fneg) = names(F1) = 
    names(opt) = names(runtime) = reg.names
  for (j in 1:N) {
    err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] = risk[[j]] =
      nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] = runtime[[j]] =
      matrix(NA,nrep,1)
  }
  filled = rep(FALSE,N)
  err.null = risk.null = sigma = rep(NA,nrep)

  # Loop through the repetitions
  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Simulation %i (of %i) ...\n",i,nrep))
      cat("  Generating data ...\n")
    }

    # Generate x, y, xval, yval
    xy.obj = sim.xy(n,p,nval,rho,s,beta.type,snr)
    risk.null[i] = diag(t(xy.obj$beta) %*% xy.obj$Sigma %*% xy.obj$beta)
    err.null[i] = risk.null[i] + xy.obj$sigma^2
    sigma[i] = xy.obj$sigma

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

        # Grab the estimated coefficients, and the predicted values on the
        # training and validation sets
        betahat = as.matrix(coef(reg.obj))
        m = ncol(betahat); nc = nrow(betahat)

        # Check for intercept
        if (nc == p+1) {
          intercept = TRUE
          betahat0 = betahat[1,]
          betahat = betahat[-1,]
        }
        else intercept = FALSE

        muhat.train = as.matrix(predict(reg.obj,xy.obj$x))
        muhat.val = as.matrix(predict(reg.obj,xy.obj$xval))

        # Populate empty matrices for our metrics, of appropriate dimension
        if (!filled[j]) {
          err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] =
            risk[[j]] = nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] =
            matrix(NA,nrep,m)
          filled[j] = TRUE
          # N.B. Filling with NAs is important, because the filled flag could
          # be false for two reasons: i) we are at the first iteration, or ii)
          # we've failed in all previous iters to run the regression method
        }

        # Record all of our metrics
        err.train[[j]][i,] = colMeans((muhat.train - xy.obj$y)^2)
        err.val[[j]][i,] = colMeans((muhat.val - xy.obj$yval)^2)
        delta = betahat - xy.obj$beta
        risk[[j]][i,] = diag(t(delta) %*% xy.obj$Sigma %*% delta)
        if (intercept) risk[[j]][i,] = risk[[j]][i,] + betahat0^2
        err.test[[j]][i,] = risk[[j]][i,] + xy.obj$sigma^2
        prop[[j]][i,] = 1 - err.test[[j]][i,] / err.null[i]
        nzs[[j]][i,] = colSums(betahat!=0)
        tpos = colSums((betahat!=0)*(xy.obj$beta!=0))
        fpos[[j]][i,] = nzs[[j]][i,]-tpos
        fneg[[j]][i,] = colSums((betahat==0)*(xy.obj$beta!=0))
        F1[[j]][i,] = 2*tpos/(2*tpos+fpos[[j]][i,]+fneg[[j]][i,])
        opt[[j]][i,] = (err.test[[j]][i,] - err.train[[j]][i,]) /
          err.train[[j]][i,]
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
      saveRDS(enlist(err.train,err.val,err.test,err.null,prop,risk,risk.null,
                     nzs,fpos,fneg,F1,opt,sigma,runtime),file=file)
    }
  }

  # Save results now (in case of an error that might occur below)
  out = enlist(err.train,err.val,err.test,err.null,prop,risk,risk.null,nzs,fpos,
               fneg,F1,opt,sigma,runtime)
  if (!is.null(file)) saveRDS(out, file)

  # Tune according to validation error, and according to test error
  out = choose.tuning.params(out)

  # Save final results
  out = c(out,list(rho=rho,s=s,beta.type=beta.type,snr=snr,call=this.call))
  class(out) = "sim"
  if (!is.null(file)) { saveRDS(out, file); invisible(out) }
  else return(out)
}

##############################

choose.tuning.params = function(obj) {
  N = length(obj$err.test) # Number of methods
  nrep = nrow(obj$err.test[[1]]) # Number of repetitions
  tun.val = matrix(NA,nrep,N)
  tun.ora = rep(NA,N)

  # Validation tuning: based on validation error
  for (i in 1:nrep) {
    for (j in 1:N) {
      tun.val[i,j] = which.min(obj$err.val[[j]][i,])
      if (length(tun.val[i,j]) == 0) tun.val[i,j] = NA
    }
  }

  # Oracle tuning: based on average test error
  for (j in 1:N) {
    tun.ora[j] = which.min(colMeans(obj$err.test[[j]], na.rm=TRUE))
  }

  return(c(obj,enlist(tun.val,tun.ora)))
}

#' Tuning and aggregation function for sim object.
#'
#' Tune and aggregate a given metric across a set of simulations, stored in a
#'   sim object (produced by \code{\link{sim.master}}).
#'
#' @param obj The sim object.
#' @param z The metric to be tuned. Should be a list of lengh N, where N is the
#'   number of methods under consideration in the simulation, and in the ith
#'   element of the list should be a matrix of dimension nrep x m, where nrep is
#'   the number of repetitions in the simulation, and m is the number of tuning
#'   parameters for the ith method.
#' @param tune Either TRUE or FALSE, indicating whether the metric should be
#'   tuned; if FALSE, then it is only aggregated (for each value of the tuning
#'   parameter). Default is TRUE.
#'
#' @return A list with components z.ave, z.std, z.med, z.mad, z.val, z.ora,
#'   z.val.ave, z.val.std, z.val.med, z.val.mad, z.ora.ave, z.ora.std,
#'   z.ora.med, z.ora.mad.  The elements z.ave, z.std, z.med, z.mad are lists of
#'   length N, where N is the number of methods in the simulation; in the ith
#'   element of z.ave is a vector of length m, the number of tuning parameters
#'   for the ith method, giving the averages of the metric across the
#'   repetitions for each tuning parameter value; similarly z.std contains
#'   standard errors, z.med contains medians, and z.mad contains median absolute
#'   deviations.  The elements z.val, z.ora are matrices of dimension nrep x N,
#'   where nrep is the number of repetitions, whose rows contain the metric
#'   after it has been tuned on each repetition by selecting the tuning
#'   parameter for each method in one of two ways: validation tuning, where the
#'   tuning parameter is selected to minimize prediction error on the validation
#'   set, and oracle tuning, where the tuning parameter is selected to minimize
#'   average test error. The elements named z.val.xxx, where xxx is one of ave,
#'   std, med, or mad, are each vectors of length N, containing the average of
#'   the metric for each method under validation tuning; similarly for the
#'   postfixes std, med, and mad; and for the elements named z.ora.xxx.
#'
#' @export tune.and.aggregate

tune.and.aggregate = function(obj, z, tune=TRUE) {
  N = length(obj$err.test) # Number of methods
  nrep = nrow(obj$err.test[[1]]) # Number of repetitions

  # Aggregate
  z.ave = z.std = z.med = z.mad = vector(mode="list",length=N)
  for (j in 1:N) {
    z[[j]] = matrix(z[[j]], nrow=nrep) # Just in case it is not a matrix
    z.ave[[j]] = colMeans(z[[j]], na.rm=TRUE)
    z.std[[j]] = apply(z[[j]], 2, sd, na.rm=TRUE) /
      sqrt(colSums(!is.na(z[[j]])))
    z.med[[j]] = apply(z[[j]], 2, median, na.rm=TRUE)
    z.mad[[j]] = apply(z[[j]], 2, mad, na.rm=TRUE) /
      sqrt(colSums(!is.na(z[[j]])))
  }
  out = enlist(z.ave,z.std,z.med,z.mad)

  # Tune and aggregate
  if (tune) {
    z.val = z.ora = vector(mode="list",length=N)
    z.val.ave = z.val.std = z.val.med = z.val.mad = rep(NA,N)
    z.ora.ave = z.ora.std = z.ora.med = z.ora.mad = rep(NA,N)

    for (j in 1:N) {
      # Validation tuning
      z.val[[j]] = z[[j]][1:nrep + (obj$tun.val[,j]-1)*nrep]
      z.val.ave[j] = mean(z.val[[j]], na.rm=TRUE)
      z.val.std[j] = sd(z.val[[j]], na.rm=TRUE) / sqrt(sum(!is.na(z.val[[j]])))
      z.val.med[j] = median(z.val[[j]], na.rm=TRUE)
      z.val.mad[j] = mad(z.val[[j]], na.rm=TRUE) / sqrt(sum(!is.na(z.val[[j]])))

      # Oracle tuning
      z.ora[[j]] = z[[j]][,obj$tun.ora[j]]
      z.ora.ave[j] = mean(z.ora[[j]], na.rm=TRUE)
      z.ora.std[j] = sd(z.ora[[j]], na.rm=TRUE) / sqrt(sum(!is.na(z.ora[[j]])))
      z.ora.med[j] = median(z.ora[[j]], na.rm=TRUE)
      z.ora.mad[j] = mad(z.ora[[j]], na.rm=TRUE) / sqrt(sum(!is.na(z.ora[[j]])))
    }

    out = c(out, enlist(z.val,z.ora,z.val.ave,z.val.std,z.val.med,z.val.mad,
                        z.ora.ave,z.ora.std,z.ora.med,z.ora.mad))
  }

  return(out)
}

##############################

#' Print function for sim object.
#'
#' Summarize and print the results of a set of simulations, stored in an object
#'   of class sim (produced by \code{\link{sim.master}}).
#'
#' @param x The sim object.
#' @param what Either "error" or "risk", indicating whether the relative test
#'     error (to the Bayes error) or relative risk (to the null risk) should be
#'     displayed. Default is "error".
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

print.sim = function(x, what=c("error","risk"), type=c("ave","med"), std=TRUE,
                     digits=3, ...) {
  what = match.arg(what)
  type = match.arg(type)

  if (!is.null(x$call)) {
    cat("\nCall:\n")
    dput(x$call)
  }

  # Construct relative test error (to Bayes) or relative risk (to null)
  N = length(x$err.test) # Number of methods
  err.rel = vector(mode="list",length=N)
  for (j in 1:N) {
    if (what=="error") err.rel[[j]] = x$err.test[[j]] / x$sigma^2
    else err.rel[[j]] = x$risk[[j]] / x$risk.null
  }
  err.obj = tune.and.aggregate(x, err.rel)
  nzs.obj = tune.and.aggregate(x, x$nzs)

  cat("\nResults for tuning parameters chosen based on validation set:\n\n")
  if (type=="ave") {
    col1 = err.obj$z.val.ave
    col2 = nzs.obj$z.val.ave
    col1.std = err.obj$z.val.std
    col2.std = err.obj$z.val.std
  }
  else {
    col1 = err.obj$z.val.med
    col2 = nzs.obj$z.val.med
    col1.std = err.obj$z.val.mad
    col2.std = err.obj$z.val.mad
  }

  tab = round(cbind(col1,col2),digits)
  tab.std = round(cbind(col1.std,col2.std), digits)
  if (std) tab = matrix(paste0(tab," (",tab.std,")"),ncol=2)
  rownames(tab) = names(x$err.test)
  colnames(tab) = c(paste("Relative",ifelse(what=="error","test error","risk")),
                    "Number of nonzeros")
  print(tab,quote=F)

  cat("\nResults for tuning parameters chosen based on test set (oracle):\n\n")
    if (type=="ave") {
    col1 = err.obj$z.ora.ave
    col2 = nzs.obj$z.ora.ave
    col1.std = err.obj$z.ora.std
    col2.std = err.obj$z.ora.std
  }
  else {
    col1 = err.obj$z.ora.med
    col2 = nzs.obj$z.ora.med
    col1.std = err.obj$z.ora.mad
    col2.std = err.obj$z.ora.mad
  }

  tab = round(cbind(col1,col2),digits)
  tab.std = round(cbind(col1.std,col2.std), digits)
  if (std) tab = matrix(paste0(tab," (",tab.std,")"),ncol=2)
  rownames(tab) = names(x$err.test)
  colnames(tab) = c(paste("Relative",ifelse(what=="error","test error","risk")),
                    "Number of nonzeros")
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
#'   is to 1:length(x$err.test), which plots all methods.
#' @param method.names the names of the methods that should be plotted. Default
#'   is NULL, in which case the names are extracted from the sim object.
#' @param what Either "error" or "risk", indicating whether the relative test
#'     error (to the Bayes error) or relative risk (to the null risk) should be
#'     displayed. Default is "error".
#' @param type Either "ave" or "med", indicating whether the average or median
#'   of the relative test error metric should be displayed. Default is "ave".
#' @param std Should standard errors be displayed (in parantheses)? When type
#'   is set to "med", the median absolute deviations are shown in place of the
#'   standard errors. Default is TRUE.
#' @param lwd,pch,main,legend graphical parameters.
#' @param make.pdf Should a pdf be produced? Default is FALSE.
#' @param fig.dir,file.name The figure directory and file name to use, only
#'   when make.pdf is TRUE. Defaults are "." and "sim". (An extension of "pdf"
#'   is always appended to the given file name.)
#' @param w,h the width and height (in inches) for the plot, used only when
#'   make.pdf is TRUE. Defaults are 6 for both.
#'
#' @export plot.sim
#' @export

plot.sim = function(x, method.nums=1:length(x$err.test), method.names=NULL,
                    what=c("error","risk"), type=c("ave","med"), std=TRUE,
                    lwd=1, pch=19, main=NULL, legend=TRUE, make.pdf=FALSE,
                    fig.dir=".", file.name="sim", w=6, h=6) {

  # Check for ggplot2 package
  if (!require("ggplot2",quietly=TRUE)) {
    stop("Package ggplot2 not installed (required here)!")
  }

  what = match.arg(what)
  type = match.arg(type)
  if (is.null(method.names)) method.names = names(x$err.test[method.nums])
  ii = method.nums

  # Construct relative test error or relative risk
  N = length(x$err.test) # Number of methods
  err.rel = vector(mode="list",length=N)
  for (j in 1:N) {
    if (what=="error") err.rel[[j]] = x$err.test[[j]] / x$sigma^2
    else err.rel[[j]] = x$risk[[j]] / x$risk.null
  }
  err.obj = tune.and.aggregate(x, err.rel)
  nzs.obj = tune.and.aggregate(x, x$nzs)

  if (type=="ave") {
    xlist = nzs.obj$z.ave
    ylist = err.obj$z.ave
    ybars = err.obj$z.std
  }
  else {
    xlist = nzs.obj$z.med
    ylist = err.obj$z.med
    ybars = err.obj$z.mad
  }

  # Produce the plot
  dat = data.frame(x=unlist(xlist[ii]),
                   y=unlist(ylist[ii]),
                   se=unlist(ybars[ii]),
                   Method=factor(rep(method.names,
                                      lapply(xlist[ii],length))))

  gp = ggplot(dat, aes(x=x,y=y,color=Method)) +
    xlab("Number of nonzero coefficients") +
    ylab(ifelse(what=="error",
                "Relative test error (to Bayes)",
                "Relative risk (to null model)")) +
    geom_line(lwd=lwd) + geom_point(pch=pch) + theme_bw()
  if (std) gp = gp + geom_errorbar(aes(ymin=y-se,ymax=y+se), width=0.2)
  if (!is.null(main)) gp = gp + ggtitle(main)
  if (!legend) gp = gp + theme(legend.pos="none")
  if (make.pdf) ggsave(sprintf("%s/%s.pdf",fig.dir,file.name),
                       height=h, width=w, device="pdf")
  else gp
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
