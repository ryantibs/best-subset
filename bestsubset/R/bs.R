#' Best subset selection.
#'
#' Compute best subset selection solutions.
#'
#' @param x Matrix of predictors, of dimension (say) n x p.
#' @param y Vector of responses, of length (say) n.
#' @param k Sparsity level, i.e., number of nonzero coefficients to allow in the
#'   subset regression model; can be a vector, in which case the best subset
#'   selection problem is solved for every value of the sparsity level. Default
#'   is 0:min(n-1,p,200) for models with intercept and 0:min(n,p,200) for models
#'   without it.
#' @param intercept Should an intercept be included in the regression model?
#'   Default is TRUE.
#' @param form Either of 1 or 2, specifying the formulation to use for best
#'   subset solution as a mixed integer quadratic program. Formulations 1 and 2
#'   correspond to equations (2.5) and (2.6) in Bertsimas, King, and Mazumder
#'   (2016), see references below. Default is 1 if n >= p, and 2 if n < p.
#' @param time.limit The maximum amount of time (in seconds) to allow Gurobi to
#'   compute the subset selection solution at each value of k. Default is 100.
#' @param params A list of parameters to pass to Gurobi. Default is list() which
#'   means no additional parameters are passed. Examples: use list(Threads=4) to
#'   specificy that Gurobi should use a maximum of 4 threads when available, or
#'   list(LogFile=f) to specify that Gurobi should print its progress to a file
#'   f. Note that if TimeLimit is an element of the params list, then its value
#'   is overriden by the last argument time.limit.
#' @param nruns The number of runs of projected gradient descent to use, where
#'   each run begins at a random initialization for the coefficients. The best
#'   estimate over all these runs (achieving the lowest criterion value) is
#'   passed to Gurobi as a warm start. Default is 50.
#' @param maxiter The maximum number of iterations for each run of projected
#'   gradient descent. Default is 1000.
#' @param tol The tolerance for the relative difference in coefficients between
#'   iterations of projected gradient descent (the algorithm terminates when the
#'   relative difference is less than the specified tolerance). Default is 1e-4.
#' @param polish Should the project gradient descent algorithm replace the
#'   estimate at each iteration by the least squares solution on the active set?
#'   Default is TRUE.
#' @param Should intermediate progress be printed out? Default is FALSE.
#'
#' @return A list with the following components:
#'   \itemize{
#'   \item beta: matrix of regression coefficients, one column per sparsity
#'     level
#'   \item status: vector of status strings returned by Gurobi's MIO solver,
#'     one element for each sparsity level
#'   \item k: vector of sparsity levels
#'   \item x, y: the passed x and y
#'   \item bx, by: the means of the columns of x, and the mean of y
#'   \item intercept: was an intercept included?
#'   }
#'
#' @details This function solves best subset selection program:
#'   \deqn{\min_\beta \|Y - X \beta\|_2^2 \;\;{\rm s.t.}\;\; \|\beta\|_0 \leq k}
#'   for a response vector \eqn{Y} and predictor matrix \eqn{X}. It follows the
#'   algorithm suggested by Bertismas, King, and Mazumder (2016) (see below for
#'   for the full reference): it first uses projected gradient descent to find
#'   an approximate solution to the above nonconvex program, and then calls
#'   Gurobi's MIO (mixed integer optimization) solver with this approximate
#'   solution as a warm start.
#'   There are two options for how to formulate best subset selection as a mixed
#'   integer quadratic program, one recommended when n >= p, and the other when
#'   n < p. They correspond to equations (2.5) and (2.6) in the paper by
#'   Bertismas, King, and Mazumder (2016), respectively.
#'
#' @author Ryan Tibshirani
#' @references This function utilizes the MIO formulation for subset selection
#'   as described in "Best subset selection via a modern optimization lens" by
#'   Dimitris Bertsimas, Angela King, and Rahul Mazumder, Annals of Statistics,
#'   44(2), 813-852, 2016. This R implementation is based on Matlab code written
#'   by Rahul Mazumder.
#' @example examples/ex.fs.R
#' @export bs

bs = function(x, y, k=0:min(nrow(x)-intercept,ncol(x),200), intercept=TRUE,
              form=ifelse(nrow(x)<ncol(x),2,1), time.limit=100, params=list(),
              nruns=50, maxiter=1000, tol=1e-4, polish=TRUE, verbose=FALSE) {

  # Check for Matrix package
  if (!require("Matrix",quietly=TRUE)) {
    stop("Package Matrix not installed (required here)!")
  }
  # Check for gurobi package
  if (!require("gurobi",quietly=TRUE)) {
    stop("Package gurobi not installed (required here)!")
  }

  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)

  # Check input data
  check.xy(x=x,y=y)

  # Save original x and y
  x0 = x
  y0 = y

  # Center and scale, etc.
  obj = standardize(x,y,intercept,FALSE)
  x = obj$x
  y = obj$y
  bx = obj$bx
  by = obj$by
  sx = obj$sx
  sy = obj$sy

  # Compute largest eigenvalue of X^T X (for step size in projected gradient)
  xtx = crossprod(x)
  if (nruns == 0) L = NULL
  else {
    if (verbose) {
      cat("0. Computing max eigenvalue of X^T X, for the")
      cat(" step size in projected gradient descent.")
    }
    L = power.method(xtx)$val
  }

  # Trim sparsity levels if we need to, and initialize some variables
  k = k[k >= 0 & k <= p]
  beta0 = NULL
  nk = length(k)
  beta = matrix(0,p,nk)
  status = rep("",nk)

  # Iterate over k vector
  for (i in 1:nk) {
    if (verbose) {
      cat(sprintf("\n%i. Solving best subset selection with k=%i.",i,k[i]))
    }

    bs.obj = bs.one.k(x,y,k[i],xtx,form=form,time.limit=time.limit,beta0=beta0,
                      L=L,nruns=nruns,maxiter=maxiter,tol=tol,polish=polish,
                      verbose=verbose,params=params)
    beta[,i] = bs.obj$beta
    status[i] = bs.obj$status
    beta0 = bs.obj$beta # Use as a warm start for the next value of k
  }
  if (verbose) cat("\n")

  # Assign column names
  colnames(beta) = as.character(k)

  out = list(beta=beta,status=status,k=k,x=x0,y=y0,bx=bx,by=by,
             intercept=intercept)
  class(out) = "bs"
  return(out)
}

# @export
bs.one.k = function(x, y, k, xtx, form=ifelse(nrow(x)<ncol(x),2,1),
                    time.limit=100, nruns=50, maxiter=1000, tol=1e-4,
                    polish=TRUE, beta0=NULL, L=NULL, verbose=FALSE,
                    params=list()) {
  n = nrow(x)
  p = ncol(x)

  # Take care of a trivial case, if needed
  if (k==0) return(list(beta=rep(0,p), status="OPTIMAL"))

  # Run the projected gradient method, gather some info from it
  best.beta = bs.proj.grad(x,y,k,nruns=nruns,maxiter=maxiter,tol=tol,
                           beta0=beta0,polish=polish,L=L,verbose=verbose)
  bigm = 2*max(abs(best.beta))
  zvec = as.numeric(best.beta != 0)

  # Set up and run the MIO solver from Gurobi. The general form is
  #   min         x^T Q x + c^T x
  #   subject to  Ax <= b
  #               l <= x <= u
  #               some x_i's binary or integral
  if (verbose) cat("\n  b. Running Gurobi's mixed integer program solver ... ")
  I = Diagonal(p,1)
  model = list()

  if (form==1) {
    rvec = c(rep(0,p), rep(1,p))
    model$A = suppressMessages(
      rbind(cbind(I,-bigm*I), cbind(-I,-bigm*I), rvec))
    model$sense = rep("<=",2*p+1)            # Ineq or eq between Ax and b?
    model$rhs = c(rep(0,2*p), k)             # The vector b
    model$ub = c(rep(bigm,p), rep(1,p))
    model$lb = c(rep(-bigm,p), rep(0,p))
    model$obj = c(-2*t(x)%*%y, rep(0,p))     # The vector c in the objective
    model$Q = bdiag(xtx, Matrix(0,p,p))
    model$vtypes = c(rep("C",p), rep("B",p)) # Variable type: cont or binary
    model$start = c(best.beta, zvec)         # Warm start from proj gradient
  }
  else {
    rvec = c(rep(0,p), rep(1,p), rep(0,n))
    model$A = suppressMessages(
      rbind(cbind(I,-bigm*I,Matrix(0,p,n)), cbind(-I,-bigm*I,Matrix(0,p,n)),
            rvec, cbind(x, Matrix(0,n,p), -Diagonal(n,1))))
    model$sense = c(rep("<=",2*p+1), rep("=",n)) # Ineq or eq between Ax and b?
    model$rhs = c(rep(0,2*p), k, rep(0,n))       # The vector b
    zeta.bd = max(colSums(apply(abs(x),1,sort,decreasing=TRUE)[1:k,,drop=F]))*bigm
    model$ub = c(rep(bigm,p), rep(1,p), rep(zeta.bd,n))
    model$lb = c(rep(-bigm,p), rep(0,p), rep(-zeta.bd,n))
    model$obj = c(-2*t(x)%*%y, rep(0,p+n))       # The vector c in the objective
    model$Q = bdiag(Matrix(0,2*p,2*p), Diagonal(n,1))
    model$vtypes = c(rep("C",p), rep("B",p), rep("C",n)) # Variable type
    model$start = c(best.beta, zvec, x%*%best.beta)      # Warm start
  }

  params$TimeLimit = time.limit
  gur.obj = quiet(gurobi(model,params))
  if (verbose) cat(sprintf("Return status: %s.", gur.obj$status))
  if (is.null(gur.obj$x)) gur.obj$x = best.beta

  return(list(beta=gur.obj$x[1:p], status=gur.obj$status))
}

bs.proj.grad = function(x, y, k, nruns=50, maxiter=1000, tol=1e-4, beta0=NULL,
                        polish=TRUE, L=NULL, verbose=FALSE) {
  n = nrow(x)
  p = ncol(x)

  # If beta0 is NULL, use thresholded least squares coefficients when p < n,
  # and thresholded marginal regression coefficients when p >= n
  if (is.null(beta0)) {
    if (p < n) beta0 = lsfit(x,y,int=FALSE)$coef
    else beta0 = t(x)%*%y/colSums(x^2)
    ids = order(abs(beta0), decreasing=TRUE)
    beta0[-ids[1:k]] = 0
  }

  # If L is NULL, use the power method to approximate the largest eigenvalue
  # of X^T X, for the step size
  if (is.null(L)) L = power.method(crossprod(x))$val

  beta.beta = beta0
  best.crit = Inf
  beta = beta0

  if (verbose) cat("\n  a. Performing projected gradient runs: ")
  for (r in 1:nruns) {
    if (verbose && (r==1 || r %% 10 ==0)) cat(sprintf("%s ... ",r))
    for (i in 1:maxiter) {
      beta.old = beta

      # Take gradient descent step
      grad = -t(x) %*% (y - x %*% beta)
      beta = beta - grad/L

      # Set to zero all but the top k
      ids = order(abs(beta), decreasing=TRUE)
      beta[-ids[1:k]] = 0

      # Perform least squares polishing, if we are asked to
      if (polish) beta[ids[1:k]] = lsfit(x[,ids[1:k]],y,int=FALSE)$coef

      # Stop if the relative difference in coefficients is small enough
      if (norm(beta - beta.old) / max(norm(beta),1) < tol) break
    }

    # Compute the criterion for the current coefficients, compare to the
    # best so far
    cur.crit = sum((y - x%*%beta)^2)
    if (cur.crit < best.crit) {
      best.crit = cur.crit
      best.beta = beta
    }

    # Start the next run off at a random spot (particular choice matches
    # Rahul's Matlab code)
    beta = beta0 + 2*runif(p)*max(abs(beta0),1)
  }

  return(best.beta)
}

##############################

quiet = function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

power.method = function(A, maxiter=100) {
  b = rnorm(ncol(A))
  for (i in 1:maxiter) {
    v = A %*% b
    b = v / norm(v)
  }
  return(list(val=norm(v)/norm(b), vec=b))
}

norm = function(v) return(sqrt(sum(v^2)))

##############################

#' Coefficient function for bs object.
#'
#' Compute coefficients at a particular sparsity level of the best subset
#'   selection model.
#'
#' @param object The bs object, as produced by the bs function.
#' @param s The sparsity level (or vector of sparsity levels) at which
#'   coefficients should be computed. If missing, then the default is use
#'   all sparsity levels of the passed bs object.
#' @param ... Other arguments (currently not used).
#'
#' @export coef.bs
#' @export

coef.bs = function(object, s, ...) {
  if (missing(s)) s = object$k
  if (any(!(s %in% object$k))) {
    stop(sprintf("s must be a subset of object$k."))
  }
  mat = matrix(rep(object$k,length(s)),nrow=length(s),byrow=TRUE)
  ind = max.col(mat==s,ties.method="first")

  beta.mat = object$beta[,ind,drop=FALSE]
  if (object$intercept) return(rbind(rep(object$by,ncol(beta.mat)),beta.mat))
  else return(beta.mat)
}

#' Predict function for bs object.
#'
#' Predict the response from a new set of predictor variables, using the
#'   coefficients from a particular step of the forward stepwise path.
#'
#' @param object The vs path object, as produced by the vs function.
#' @param newx Matrix of new predictor variables at which predictions should
#'   be made; if missing, the original (training) predictors are used.
#' @param s The sparsity level (or vector of sparsity levels) at which
#'   coefficients should be computed. If missing, then the default is use
#'   all sparsity levels of the passed bs object.
#' @param ... Other arguments (currently not used).
#'
#' @export predict.bs
#' @export

predict.bs = function(object, newx, s, ...) {
  beta = coef.bs(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))

  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}
