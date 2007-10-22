xyw.coords <- function(x, y = NULL, w = NULL, v = NULL, ...) {
  xy <- xy.coords(x, y)
  x <- xin <- xy$x
  y <- yin <- xy$y
  m <- length(x)
  if(is.null(w)) w <- rep(1, m)
  if(length(w) == 1) w <- rep(w, m)
  if(length(w) != m)
    stop("length of 'w' (weights) needs to be ", m, ", not ", length(w))
  if(is.null(v)) v <- rep(0, m)
  if(length(v) == 1) v <- rep(v, m)
  if(length(v) != m)
    stop("length of 'w' (weights) needs to be ", m, ", not ", length(v))
  x <- signif(x, 6)
  if(any(duplicated(x))) {
    ux <- unique(sort(x))
    ox <- match(x, ux)
    y <- y * w
    y <- tapply(y, ox, sum)
    w <- tapply(w, ox, sum)
    v <- tapply(v, ox, sum)
    v <- ifelse(v > 0, 1, ifelse(v < 0, -1, 0))
    y <- y/w
    x <- tapply(x, ox, mean)
  }
  m <- length(x)
  if(m < 4) stop("number of observations must be at least 4")
  o <- order(x)
  list(x = x[o], y = y[o], w = w[o], v = v[o], xin = xin, yin = yin)
}
