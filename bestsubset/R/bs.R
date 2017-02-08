#' Best subset selection.
#'
#' Compute best subset selection solutions.
#'
#' @param x Matrix of predictors, of dimension (say) n x p.
#' @param y Vector of responses, of length (say) n.
#' @param k Sparsity level, i.e., number of nonzero coefficients to allow in the
#'   subset regression model; can be a vector, in which case the best subset
#'   selection problem is solved for every value of the sparsity level. Default
#'   is 1:min(n,p).
#' @param intercept Should an intercept be included in the regression model? 
#'   Default is TRUE.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
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
#'   \deqn{\min_\beta \|Y - X \beta\|_2^2 \;\;\text{subject to}\;\;
#'     \|\beta\|_0 \leq k}
#'   for a response vector \eqn{Y} and predictor matrix \eqn{X}. It uses
#'   projected gradient descent to find an approximate solution to the
#'   above nonconvex program, and then calls Gurobi's MIO (mixed integer
#'   optimization) solver with this approximate solution as a warm start.
#'   See references below for the paper by Bertsimas, King, and Mazumder
#'   (2016), that describes this algorithm.
#' 
#' @author Ryan Tibshirani 
#' @references This function utilizes the MIO formulation for subset selection
#'   as described in "Best subset selection via a modern optimization lens" by
#'   Dimitris Bertsimas, Angela King, and Rahul Mazumder, Annals of Statistics,
#'   44(2), 813-852, 2016. This R implementation is based on Matlab code written
#'   by Rahul Mazumder.
#' @example examples/ex.bs.R
#' @export bs

bs = function(x, y, k=1:min(nrow(x),ncol(x)), intercept=TRUE, time.limit=100,
              nruns=50, maxiter=1000, tol=1e-4, polish=TRUE, verbose=FALSE) {
  
  # Check for Matrix package
  if (!require("Matrix",quietly=TRUE)) {
    stop("Package Matrix not installed (required here)!")
  }
  # Check for gurobi package
  if (!require("gurobi",quietly=TRUE)) {
    stop("Package Matrix not installed (required here)!")
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
  
  beta0 = NULL
  nk = length(k)
  beta = matrix(0,p,nk)
  status = rep("",nk)
  
  # Iterate over k vector
  for (i in 1:nk) {
    if (verbose) {
      cat(sprintf("\n%i. Solving best subset selection with k=%i.",i,k[i]))
    }
    
    bs.obj = bs.one.k(x,y,k[i],xtx,time.limit=time.limit,beta0=beta0,L=L,
                      nruns=nruns,maxiter=maxiter,tol=tol,polish=polish,
                      verbose=verbose)
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
  
bs.one.k = function(x, y, k, xtx, time.limit=100, nruns=50, maxiter=1000,
                    tol=1e-4, polish=TRUE, beta0=NULL, L=NULL, verbose=FALSE) {       
  n = nrow(x)
  p = ncol(x)
  
  # Run the projected gradient method, gather some info from it
  best.beta = bs.proj.grad(x,y,k,nruns=nruns,maxiter=maxiter,tol=tol,
                           polish=polish,L=L,verbose=verbose)
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
  model$A = rbind(cbind(I, -bigm*I), cbind(-I, -bigm*I), c(rep(0,p),rep(1,p)))
  model$sense = c(rep("<=",2*p),"=")       # The constraints between Ax and b
  model$rhs = c(rep(0,2*p),k)              # The vector b
  model$ub = c(rep(bigm,p), rep(1,p)) 
  model$lb = c(rep(-bigm,p), rep(0,p))
  model$obj = c(-2*t(x)%*%y, rep(0,p))     # The vector c in the objective
  model$Q = bdiag(xtx, Matrix(0,p,p))
  model$vtypes = c(rep("C",p), rep("B",p)) # Variable type: continuous or binary
  model$start = c(best.beta, zvec)         # Warm start, best proj gradient run
  
  params = list()
  if (!is.null(time.limit)) params$TimeLimit = time.limit
  gur.obj = quiet(gurobi(model,params))

  if (verbose) cat(sprintf("Return status: %s.", gur.obj$status))
  
  return(list(beta=gur.obj$x[1:p], status=gur.obj$status))
}

bs.proj.grad = function(x, y, k, nruns=50, maxiter=1000, tol=1e-4, polish=TRUE,
                        beta0=NULL, L=NULL, verbose=FALSE) {
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
    if (cur.crit < best.crit) best.beta = beta
    
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
#' @export 

coef.bs = function(object, s, ...) {
  if (missing(s)) s = object$k
  else if (any(!(s %in% object$k))) {
    stop(sprintf("s must be a subset of object$k."))
  }
  return(object$beta[,s])           
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
#' @export 

predict.bs = function(object, newx, s, ...) {
  beta = coef.bs(object,s)
  if (missing(newx)) newx = object$x
  else {
    newx = matrix(newx,ncol=ncol(object$x))
    newx = scale(newx,object$bx,FALSE)
  }
  return(newx %*% beta + object$by)
}

##############################

# Interpolation function to get coefficients

coef.interpolate = function(beta, s, knots, decreasing=TRUE) {
  # Sort the s values 
  o = order(s,decreasing=decreasing)
  s = s[o]
        
  k = length(s)
  mat = matrix(rep(knots,each=k),nrow=k)
  if (decreasing) b = s >= mat
  else b = s <= mat
  blo = max.col(b,ties.method="first")
  bhi = pmax(blo-1,1)

  i = bhi==blo
  p = numeric(k)
  p[i] = 0
  p[!i] = ((s-knots[blo])/(knots[bhi]-knots[blo]))[!i]
  
  beta = t((1-p)*t(beta[,blo,drop=FALSE]) + p*t(beta[,bhi,drop=FALSE]))
  colnames(beta) = as.character(round(s,3))
  rownames(beta) = NULL

  # Return in original order
  o = order(o)
  return(beta[,o,drop=FALSE])
}
