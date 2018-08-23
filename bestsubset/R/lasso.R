#' Lasso and friends.
#'
#' Compute the lasso, ridge regression, or elastic net solutions in regression.
#'
#' @param nrelax The number of interpolations to produce between the lasso
#'   (or ridge, or elastic net) solution, at each value of lambda, and the
#'   least squares coefficients on the corresponding active set. (The number
#'   of interpolations counts the endpoints inclusively, i.e., the lasso and
#'   least squares solutions.) Default is 1, which means that only the lasso
#'   solution is considered (no strict relaxations), at each value of lambda.
#'
#' @description This is just a simple wrapper function around the
#'   \code{\link{glmnet}} function in the R package of the same name. Its
#'   purpose is twofold: (i) to provide functionality where the associated coef
#'   and predict methods always produce coefficients and predictions at exactly
#'   nlambda values, by default (the \code{\link{glmnet}} function may produce a
#'   path with less than nlambda lambda values, depending on the data); and (ii)
#'   to provide relaxed versions of the lasso (or ridge regression, or elastic
#'   net) solutions, defined by interpolating in between each solution and the
#'   least squares coefficients on the corresponding active set. The number of
#'   interpolations is governed by the nrelax argument; all other arguments are
#'   the same as in \code{\link{glmnet}}.
#'
#' @author Trevor Hastie, Robert Tibshirani, Ryan Tibshirani
#' @export lasso

lasso = function(x, y, alpha=1, nrelax=1, nlambda=50,
                 lambda.min.ratio=ifelse(nrow(x)<ncol(x),0.01,0.0001),
                 lambda=NULL, intercept=TRUE, standardize=TRUE) {
  
  # Check for glmnet package
  if (!require("glmnet",quietly=TRUE)) {
    stop("Package glmnet not installed (required here)!")
  }

  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)

  # Set dfmax manually
  dfmax = p
  if (nrelax > 1 && n < (p-intercept)) dfmax = n-intercept
  # TH: this sometimes still returns too many!
  
  # Reset nlambda if a specific lambda sequence is passed
  if (!is.null(lambda)) nlambda = length(lambda)

  # Run glmnet
  obj = glmnet(x, y, alpha=alpha, nlambda=nlambda, dfmax=dfmax,
               lambda.min.ratio=lambda.min.ratio, lambda=lambda,
               intercept=intercept, standardize=standardize)
  
  # Append a few things to the returned object
  obj$nrelax = nrelax
  obj$nlambda = nlambda
  obj$intercept = intercept
  obj$x = x; obj$y = y
  class(obj) = "lasso"
  return(obj)
}

#' Coef function for lasso object.
#' @export coef.lasso
#' @export

coef.lasso = function(object, s=NULL, gamma=NULL) {
  beta.lasso = coef.lasso.from.glmnet(object,s)
  if (object$nrelax == 1 && is.null(gamma)) {
    if (object$intercept) return(beta.lasso)
    else return(beta.lasso[-1,])
  }
  if (is.null(gamma)) gamma = seq(1,0,length=object$nrelax)
  
  beta.ls = coef.ls(beta.lasso,object$x,object$y)
  beta.left = matrix(apply(beta.lasso,2,function(b){b%o%gamma}),
                     nrow=nrow(beta.lasso))
  beta.right = matrix(apply(beta.ls,2,function(b){b%o%(1-gamma)}),
                      nrow=nrow(beta.lasso))
  beta.mat = beta.left + beta.right
  
  if (object$intercept) return(beta.mat)
  else return(beta.mat[-1,])
}

coef.lasso.from.glmnet = function(object, s=NULL) {
  class(object) = "glmnet"
  if (length(object$lambda)==object$nlambda) {
    return(glmnet::coef.glmnet(object,s=s))
  }
  else {
    min.lam = min(object$lambda)
    max.lam = max(object$lambda)
    svec = exp(seq(log(max.lam),log(min.lam),length=object$nlambda))
    return(glmnet::coef.glmnet(object,s=svec))
    ## RJT TODO: should we used exact=TRUE above? Requires additional
    ## arguments to match the initial call to glmnet(), kind of clunky
    ## TH: use glmnet.control(fdev=0) at beginning of session
    ## Still needed though for cases when df exceeds p (can happen with
    ## glmnet, and bad for relaxed lasso)
  }
}

coef.ls = function(beta, x, y) {
  n = nrow(x); p = ncol(x)
  apply(beta, 2, function(b) {
    act.set = which(b[-1] != 0)
    intercept = b[1]!=0
    if (length(act.set)==0) return(c(b[1],rep(0,p)))
    if (length(act.set)>(n-intercept)) {
      # Take any n-intercept elements (which ones dont matter)
      act.set = act.set[seq(n-intercept)]
    }
    b.new = rep(0,p+1)
    if (intercept) b.new[c(1,1+act.set)] = lsfit(x[,act.set],y)$coef
    else b.new[1+act.set] = lsfit(x[,act.set],y,int=FALSE)$coef
    return(b.new)
  })
}

#' Predict function for lasso object.
#' @export predict.lasso
#' @export

predict.lasso = function(object, newx, s=NULL) {
  if (missing(newx)) newx = object$x
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% coef.lasso(object,s))
}

