#' Forward stepwise regression.
#'
#' Compute the forward stepwise regression path.
#'
#' @param x Matrix of predictors, of dimension (say) n x p.
#' @param y Vector of responses, of length (say) n.
#' @param maxsteps Maximum number of steps of the forward stepwise path to
#'   compute. Default is min(n-1,p,2000) for models with intercept and
#'   min(n,p,2000) for models without it
#' @param intercept,normalize Should an intercept be included in the regression
#'   model? Should the predictors be normalized before computing the path?
#'   Default is TRUE for both.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#'
#' @return A list with the following components:
#'   \itemize{
#'   \item action: vector giving the index of the variable added at each step
#'   \item df: vector that gives the (naive) degrees of freedom of the model
#'   at each step, i.e., the number of active predictor variables (+ 1 if there
#'   is an intercept included)
#'   \item beta: matrix of regression coefficients for each step along the path,
#'     one column per step
#'   \item completepath: a boolean indicating whether the forward stepwise path
#'   was run to completion (as opposed to being stopped early because the max
#'   number of steps was achieved)
#'   \item bls: if the complete path was computed, this is a vector that gives
#'   the least squares coefficients of the full regression model
#'   \item x, y: the passed x and y
#'   \item bx, by: the means of the columns of x, and the mean of y
#'   \item intercept, normalize: the passed values for intercept and normalize
#'   }
#'
#' @details This function implements forward stepwise regression, adding the
#'   predictor at each step that maximizes the absolute correlation between the
#'   predictors---once orthogonalized with respect to the current model---and
#'   the residual. This entry criterion is standard, and is equivalent to
#'   choosing the variable that achieves the biggest drop in RSS at each step;
#'   it is used, e.g., by the \code{step} function in R. Note that, for example,
#'   the \code{lars} package implements a stepwise option (with type="step"),
#'   but uses a (mildly) different entry criterion, based on maximal absolute
#'   correlation between the original (non-orthogonalized) predictors and the
#'   residual.
#'
#' @author Ryan Tibshirani
#' @example examples/ex.fs.R
#' @export fs

fs = function(x, y, maxsteps=min(nrow(x)-intercept,ncol(x),2000),
              intercept=TRUE, normalize=TRUE, verbose=FALSE) {

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
  obj = standardize(x,y,intercept,normalize)
  x = obj$x
  y = obj$y
  bx = obj$bx
  by = obj$by
  sx = obj$sx

  #####
  # Find the first variable to enter and its sign
  z = scale(x,center=F,scale=sqrt(colSums(x^2)))
  u = t(z) %*% y
  j.hit = which.max(abs(u))   # Hitting coordinate
  sign.hit = Sign(u[j.hit])   # Hitting sign

  if (verbose) {
    cat(sprintf("1. Added variable %i, |A|=%i...",j.hit,1))
  }

  # Now iterate to find the sequence of FS estimates

  # Things to keep track of, and return at the end
  buf = min(maxsteps+1,500)
  action = numeric(buf)      # Actions taken
  df = numeric(buf)          # Degrees of freedom
  beta = matrix(0,p,buf)     # FS estimates

  # Record action, df, solution (df and solution are here
  # correspond to step 0; always a step behind)
  action[1] = j.hit
  df[1] = 0
  beta[,1] = 0

  # Other things to keep track of, but not return
  r = 1                       # Size of active set
  A = j.hit                   # Active set
  I = Seq(1,p)[-j.hit]        # Inactive set
  sign = sign.hit             # Active signs
  X1 = x[,j.hit,drop=FALSE]   # Matrix X[,A]
  X2 = x[,-j.hit,drop=FALSE]  # Matrix X[,I]
  k = 2                       # Step counter

  # Compute a skinny QR decomposition of X1
  qr.obj = qr(X1)
  Q = qr.Q(qr.obj,complete=TRUE)
  Q1 = Q[,1,drop=FALSE];
  Q2 = Q[,-1,drop=FALSE]
  R = qr.R(qr.obj)

  # Throughout the algorithm, we will maintain
  # the decomposition X1 = Q1*R. Dimensions:
  # X1: n x r
  # Q1: n x r
  # Q2: n x (n-r)
  # R:  r x r

  while (k <= maxsteps) {
    ##########
    # Check if we've reached the end of the buffer
    if (k > length(action)) {
      buf = length(action)
      action = c(action,numeric(buf))
      df = c(df,numeric(buf))
      beta = cbind(beta,matrix(0,p,buf))
    }

    # Key quantities for the next entry
    a = backsolve(R,t(Q1) %*% y)
    b = backsolve(R,t(Q1) %*% X2)
    X2.resid = X2 - X1 %*% b
    z = scale(X2.resid,center=F,scale=sqrt(colSums(X2.resid^2)))
    u = as.numeric(t(z) %*% y)

    # Otherwise find the next hitting time
    sign.u = Sign(u)
    abs.u = sign.u * u
    j.hit = which.max(abs.u)
    sign.hit = sign.u[j.hit]

    # Record action, df, solution
    action[k] = I[j.hit]
    df[k] = r
    beta[A,k] = a

    # Update rest of the variables
    r = r+1
    A = c(A,I[j.hit])
    I = I[-j.hit]
    sign = c(sign,sign.hit)
    X1 = cbind(X1,X2[,j.hit])
    X2 = X2[,-j.hit,drop=FALSE]

    # Update the QR decomposition
    updated.qr = updateQR(Q1,Q2,R,X1[,r])
    Q1 = updated.qr$Q1
    Q2 = updated.qr$Q2
    R = updated.qr$R

    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...",k,A[r],r))
    }

    # Update counter
    k = k+1
  }

  # Record df and solution at last step
  df[k] = k-1
  beta[A,k] = backsolve(R,t(Q1) %*% y)

  # Trim
  action = action[Seq(1,k-1)]
  df = df[Seq(1,k)]
  beta = beta[,Seq(1,k),drop=FALSE]

  # If we stopped short of the complete path, then note this
  if (k-1 < min(n-intercept,p)) {
    completepath = FALSE
    bls = NULL
  }

  # Else we computed the complete path, so record LS solution
  else {
    completepath = TRUE
    bls = beta[,k]
  }

  if (verbose) cat("\n")

  # Adjust for the effect of centering and scaling
  if (intercept) df = df+1
  if (normalize) beta = beta/sx
  if (normalize && completepath) bls = bls/sx

  # Assign column names
  colnames(beta) = as.character(Seq(0,k-1))

  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  class(out) = "fs"
  return(out)
}

##############################

# Downdate the QR factorization, after a column has
# been deleted. Here Q1 is m x n, Q2 is m x k, and
# R is n x n.
#' @useDynLib bestsubset downdate1

downdateQR = function(Q1,Q2,R,col) {
  m = nrow(Q1)
  n = ncol(Q1)

  a = .C("downdate1",
    Q1=as.double(Q1),
    R=as.double(R),
    col=as.integer(col-1),
    m=as.integer(m),
    n=as.integer(n),
    dup=FALSE,
    PACKAGE="bestsubset")

  Q1 = matrix(a$Q1,nrow=m)
  R = matrix(a$R,nrow=n)

  # Re-structure: add a column to Q2, delete one from
  # Q1, and trim R
  Q2 = cbind(Q2,Q1[,n])
  Q1 = Q1[,-n,drop=FALSE]
  R = R[-n,-col,drop=FALSE]

  return(list(Q1=Q1,Q2=Q2,R=R))
}

# Update the QR factorization, after a column has been
# added. Here Q1 is m x n, Q2 is m x k, and R is n x n.
#' @useDynLib bestsubset update1

updateQR = function(Q1,Q2,R,col) {
  m = nrow(Q1)
  n = ncol(Q1)
  k = ncol(Q2)

  a = .C("update1",
    Q2=as.double(Q2),
    w=as.double(t(Q2) %*% col),
    m=as.integer(m),
    k=as.integer(k),
    dup=FALSE,
    PACKAGE="bestsubset")

  Q2 = matrix(a$Q2,nrow=m)
  w = c(t(Q1) %*% col,a$w)

  # Re-structure: delete a column from Q2, add one to
  # Q1, and expand R
  Q1 = cbind(Q1,Q2[,1])
  Q2 = Q2[,-1,drop=FALSE]
  R = rbind(R,rep(0,n))
  R = cbind(R,w[Seq(1,n+1)])

  return(list(Q1=Q1,Q2=Q2,R=R))
}

##############################

#' Coefficient function for fs object.
#'
#' Compute coefficients at a particular step of the forward stepwise path.
#'
#' @param object The fs object, as produced by the fs function.
#' @param s The step (or vector of steps) of the path at which coefficients
#'   should be computed. Can be fractional, in which case interpolation is
#'   performed. If missing, then the default is use all steps of the passed
#'   fs object.
#' @param ... Other arguments (currently not used).
#'
#' @details Note that at s = 1, there is one nonzero coefficient, at
#'   s = 2, there are two nonzero coefficients, etc. (This differs from the
#'   parametrization used in the \code{coef.fs} function in the R package
#'   \code{selectiveInference}, as the latter function delivers s-1 nonzero
#'   coefficients at step s, and was written to be consistent with the
#'   natural parametrization for the least angle regression path.)
#'
#' @export coef.fs
#' @export

coef.fs = function(object, s, ...) {
  beta = object$beta
  if (object$completepath) beta = cbind(beta,object$bls)
  k = ncol(beta)-1
  if (missing(s)) s = 0:k
  else if (min(s)<0 || max(s)>k) stop(sprintf("s must be between 0 and %i",k))
  knots = 0:k
  decreasing = FALSE

  beta.mat = coef.interpolate(beta,s,knots,decreasing)
  if (object$intercept) return(rbind(rep(object$by,ncol(beta.mat)),beta.mat))
  else return(beta.mat)
}

#' Predict function for fs object.
#'
#' Predict the response from a new set of predictor variables, using the
#'   coefficients from a particular step of the forward stepwise path.
#'
#' @param object The fs path object, as produced by the fs function.
#' @param newx Matrix of new predictor variables at which predictions should
#'   be made; if missing, the original (training) predictors are used.
#' @param s The step (or vector of steps) of the path at which coefficients
#'   should be computed. Can be fractional, in which case interpolation is
#'   performed. If missing, then the default is use all steps of the passed
#'   fs object.
#' @param ... Other arguments (currently not used).
#'
#' @details Note that at s = 1, there is one nonzero coefficient, at
#'   s = 2, there are two nonzero coefficients, etc. (This differs from the
#'   parametrization used in the \code{coef.fs} function in the R package
#'   \code{selectiveInference}, as the latter function delivers s-1 nonzero
#'   coefficients at step s, and was written to be consistent with the
#'   natural parametrization for the least angle regression path.)
#'
#' @export predict.fs
#' @export

predict.fs = function(object, newx, s, ...) {
  beta = coef.fs(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
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
