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
#'   variables; specifically, predictors i and j will have pairwise correlation
#'   rho^abs(i-j). Default is 0.
#' @param s number of nonzero coefficients in the underlying regression model.
#'   Default is 5. (Ignored if beta.type is 4, in which case the number of
#'   nonzero coefficients is 6; and if beta.type is 5, it is interpreted as a
#'   the number of strongly nonzero coefficients in a weak sparsity model.)
#' @param beta.type Integer taking values in between 1 and 5, used to specify
#'   the pattern of nonzero coefficients in the underlying regression model; see
#'   details below. Default is 1.
#' @param snr Desired signal-to-noise ratio; the variance of the errors will be
#'   set accordingly. Default is 1.
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
    beta[1:s] = seq(10,0.5,length=10)
  } else if (beta.type==4) {
    beta[1:6] = c(-10,-6,-2,2,6,10)
  } else {
    beta[1:s] = 1
    beta[(s+1):p] = 0.8^(1:(p-s))
  }

  # Set snr based on sample variance on large test set 
  mutest = xtest %*% beta
  vmu = var(mutest) 
  sigma = sqrt(vmu/snr)

  # Generate responses
  y = x %*% beta + rnorm(n)*sigma
  yval = xval %*% beta + rnorm(nval)*sigma
  ytest = mutest + rnorm(ntest)*sigma
 
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
#' @param file,file.rep Name of a file to which simulation results will be saved
#'   (using saveRDS), and a number of repetitions after which intermediate
#'   results will be saved. Setting file to NULL is interpreted to mean that no
#'   simulations results should be saved; setting file.rep to 0 is interpreted
#'   to mean that simulations results should be saved at the very end, i.e., no
#'   intermediate saving. Defaults are NULL and 5, respectively.
#' @param rho,s,beta.type,snr Arguments to pass to \code{\link{sim.xy}}; see the
#'   latter's help file for details.
#'
#' @return A list with components err.train, err.val, err.test, risk, nzs, opt.
#'     for the training error, validation error, test error, risk, number of
#'     selected nonzero coefficients, and optimism (difference in test and
#'     training errors). These are each lists of length N, where N is the number
#'     of regression methods under consideration (the length of reg.funs). The
#'     ith element of each list is then a matrix of dimension nrep x m, where m 
#'     the number of tuning parameters inherent to the ith method. The returned
#'     components err.train.ave, err.val.ave, err.test.ave, risk.ave, nzs.ave,
#'     opt.ave return the averages of the training error, validation error, etc.
#'     over the nrep repetitions. Similarly for the components with postfixes
#'     .std, .med, .mad, which return the standard deviation, median, and median
#'     absolute deviation, respectively. The returned components with postfixes
#'     .tun.val and .tun.oracle are matrices of dimension N x nrep, which return
#'     the training error, validation error, etc. when the tuning parameter for
#'     each regression method in each repetition is chosen by validation tuning
#'     (best validation error) or by oracle tuning (best average risk across all
#'     the repetitions). 
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
  nms = names(reg.funs)
  if (is.null(nms)) nms = paste("Method",1:N)
  
  err.train = err.val = err.test = vector(mode="list",length=N)
  risk = nzs = opt = tim = vector(mode="list",length=N)
  for (j in 1:N) {
    err.train[[j]] = err.val[[j]] = err.test[[j]] = matrix(NA,nrep,1)
    risk[[j]] = nzs[[j]] = opt[[j]] = tim[[j]] = matrix(NA,nrep,1)
  }
  filled = rep(FALSE,N)

  # Loop through the repetitions
  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Simulation %i (of %i) ...\n",i,nrep))
      cat("\tGenerating data ...\n")
    }
    
    # Generate x, y, xval, yval, xtest, ytest
    xy.obj = sim.xy(n,p,nval,ntest,rho,s,beta.type,snr)

    # Loop through the regression methods
    for (j in 1:N) {
      if (verbose) {
        cat(sprintf("\tApplying regression method %i (of %i) ...\n",
                    j,N))
      }
      
      tryCatch({
        # Apply the regression method in hand
        tim[[j]][i] = system.time({
          reg.obj = reg.funs[[j]](xy.obj$x,xy.obj$y)
        })[1]

        # Grab the coefficients, and predicted values on the training,
        # validation, and testing sets
        beta = coef(reg.obj)
        muhat.train = predict(reg.obj,xy.obj$x)
        muhat.val = predict(reg.obj,xy.obj$xval)
        muhat.test = predict(reg.obj,xy.obj$xtest)
        
        # Populate empty matrices for our metrics, of appropriate dimension
        if (!filled[j]) {
          err.train[[j]] = err.val[[j]] = err.test[[j]] = 
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
        risk[[j]][i,] = colMeans((muhat.test - xy.obj$mutest)^2)
        nzs[[j]][i,] = colMeans(beta != 0)
        opt[[j]][i,] = err.test[[j]][i,] - err.train[[j]][i,]
      }, error = function(err) {
        if (verbose) {
          cat(paste("\t\tOops! Something went wrong, see error message",
                    "below; recording all metrics here as NAs ...\n"))
          cat("\t\t***** Error message *****\n")
          cat(sprintf("\t\t%s\n",err$message))
          cat("\t\t*** End error message ***\n")
        }
        # N.B. No need to do anything, the metrics are already filled with NAs
      })
    }

    # Save intermediate results?
    if (!is.null(file) && file.rep > 0 && i %% file.rep == 0) {
      saveRDS(enlist(err.train,err.val,err.test,risk,nzs,opt),file)
    }
  }

  # Save results now (in case of error below)
  res = enlist(err.train,err.val,err.test,risk,nzs,opt)
  if (!is.null(file)) saveRDS(res, file)
  
  # Aggregate our metrics over the simulations
  stats = compute.stats(res)

  # Tune according to validation error, and according to oracle
  tuned = tune.methods(res, stats)
  
  # Save final results
  out = c(res, stats, tuned)
  class(out) = "sim"
  if (!is.null(file)) saveRDS(out, file)
  return(out)
}

compute.stats = function(res) {
  m = length(res) # Number of metrics
  N = length(res$err.val) # Number of methods
  res.ave = res.std = vector(mode="list", length=m)
  res.med = res.mad = vector(mode="list", length=m)
  names(res.ave) = paste0(names(res), ".ave")
  names(res.std) = paste0(names(res), ".std")
  names(res.med) = paste0(names(res), ".med")
  names(res.mad) = paste0(names(res), ".mad")
  
  for (i in 1:m) {
    res.ave[[i]] = res.std[[i]] = res.med[[i]] = res.mad[[i]] =
      vector(mode="list", length=N)
    names(res.ave[[i]]) = names(res.std[[i]]) = names(res.med[[i]]) =
      names(res.mad[[i]]) = names(res[[i]])
    
    for (j in 1:N) {
      res.ave[[i]][[j]] = rowMeans(res[[i]][[j]], na.rm=TRUE)
      res.std[[i]][[j]] = apply(res[[i]][[j]], 2, sd, na.rm=TRUE) /
        sqrt(rowSums(!is.na(res[[i]][[j]])))
      res.med[[i]][[j]] = apply(res[[i]][[j]], 2, median, na.rm=TRUE)
      res.mad[[i]][[j]] = apply(res[[i]][[j]], 2, mad, na.rm=TRUE) /
        sqrt(rowSums(!is.na(res[[i]][[j]])))
    }
  }
  
  return(c(res.ave,res.std,res.med,res.mad))
}

tune.methods = function(res, stats) {
  # Tune first by min validation error
  N = length(res$err.val) # Number of methods
  nrep = nrow(res$err.val[[1]]) # Number of repetitions
  err.test.tun.val = nzs.tun.val = matrix(NA,nrep,N)
  
  for (i in 1:nrep) {
    for (j in 1:N) {
      k.val = which.min(res$err.val[[j]][i,])
      err.test.tun.val[i,j] = res$err.test[[j]][i,k.val]
      nzs.tun.val[i,j] = res$nzs[[j]][i,k.val]
    }
  }

  # Tune second by min oracle error (average risk)
  err.test.tun.oracle = nzs.tun.oracle = matrix(NA,nrep,N)
  
  for (j in 1:N) {
    k.oracle = which.min(stats$risk.ave[[j]])
    for (i in 1:nrep) {
      err.test.tun.oracle[i,j] = res$err.test[[j]][i,k.oracle]
      nzs.tun.oracle[i,j] = res$nzs[[j]][i,k.oracle]
    }
  }

  return(enlist(err.test.tun.val,nzs.tun.val,
                err.test.tun.oracle,nzs.tun.oracle))
}
