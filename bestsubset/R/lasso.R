#' Lasso and friends.
#'
#' Compute the lasso, ridge regression, or elastic net solutions in regression.
#'
#' @description This is just a simple wrapper function around the
#'   \code{\link{glmnet}} function in the R package of the same name. Its 
#'   purpose is to provide a version where the associated coef and predict
#'   methods always produce coefficients and predictions at exactly nlambda
#'   values, by default. (The \code{\link{glmnet}} function may produce a
#'   path with less than nlambda lambda values, depending on the data.)
#' 
#' @export lasso

lasso = function(x, y, alpha=1, nrelax=10, nlambda=50,
                 lambda.min.ratio=ifelse(nrow(x)<ncol(x),0.01,0.0001),
                 lambda=NULL, intercept=TRUE, standardize=TRUE) {

  # Check for glmnet package
  if (!require("glmnet",quietly=TRUE)) {
    stop("Package glmnet not installed (required here)!")
  }
  
  obj = glmnet(x,y,alpha=alpha,nlambda=nlambda,
               lambda.min.ratio=lambda.min.ratio,
               lambda=lambda,intercept=intercept,
               standardize=standardize)
  obj$nlambda = nlambda
  class(obj) = "lasso"
  return(obj)
}

#' Coef function for lasso object.
#' @export 

coef.lasso = function(object, s=NULL) {
  predict.lasso(object,s=s,type="coefficients")
}


#' Predict function for lasso object.
#' @export 

predict.lasso = function(object, newx, s=NULL, type="link") {
  if (length(object$lambda)==object$nlambda) {
    return(glmnet::predict.glmnet(object,newx,s=s,type=type))
  }
  else {
    min.lam = min(object$lambda)
    max.lam = max(object$lambda)
    svec = log(seq(exp(max.lam),exp(min.lam),length=object$nlambda))
    return(glmnet::predict.glmnet(object,newx,s=svec,type=type))
  }
}

