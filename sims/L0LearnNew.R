# New coef and predict functions for L0Learn. Just allow them to take no lambda
# or gamma arguments, and return coefficients/predictions over all tunings
coef.L0LearnNew = function(object) {
  class(object) = "L0Learn"
  mat = coef(object,lambda=object$lambda[[1]],gamma=object$gamma[1])
  if (length(object$gamma) > 1) {
    for (i in 2:length(object$gamma)) {
      mat = cbind(mat, coef(object, lambda=object$lambda[[i]],
                            gamma=object$gamma[i]))
    }
  }
  return(mat)
}

predict.L0LearnNew = function(object,newx) {
  class(object) = "L0Learn"
  mat = predict(object,newx,lambda=object$lambda[[1]],gamma=object$gamma[1])
  if (length(object$gamma) > 1) {
    for (i in 2:length(object$gamma)) {
      mat = cbind(mat, predict(object, newx,lambda=object$lambda[[i]],
                               gamma=object$gamma[i]))
    }
  }
  return(mat)
}
