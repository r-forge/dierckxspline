predict.dierckx <- function(object, newx = NULL, ...) {
  if(is.null(newx)) return(fitted(object))
  if(!is.numeric(newx) || length(newx) == 0)
    stop("Invalid value for argument 'newx'")
  if(!is.null(object$periodic) && object$periodic) {
    period <- range(object$x)
    xb <- period[1]
    xe <- period[2]
    if(any(newx < xb | newx > xe)) {
      rn <- range(newx) %/% diff(period)
      br <- if(rn[1] == rn[2]) {
        c(rn[1], rn[1] + diff(period))
      } else {
        seq(rn[1], rn[2], diff(period))
      }
      ct <- cut(newx, br, include.lowest = TRUE)
      newx <- newx - br[ct]
    }
  }
  o <- order(newx)
  m <- length(newx)
  sp <- .Fortran("splev",
                 knots = as.single(object$knots),
                 n = as.integer(object$n),
                 coef = as.single(object$coef),
                 k = as.integer(object$k),
                 x = as.single(newx[o]),
                 sp = single(m),
                 m = as.integer(m),
                 ier = as.integer(object$ier))$sp
  sp <- as.numeric(sp)
  sp <- sp[order(o)]
  sp
}
