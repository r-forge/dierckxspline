coef.dierckx <- function(object, ...) {
  object$coef[abs(object$coef) > 0]
}
