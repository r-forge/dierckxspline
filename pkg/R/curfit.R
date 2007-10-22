curfit <- function(x, ...) {
  UseMethod("curfit")
}

percur <- function(x, ...) {
  curfit(x, periodic = TRUE, ...)
}

curfit.formula <- function(x, data, ..., weights, subset) {
#  cl <- match.call()
#  mf <- match.call(expand.dots = FALSE)
#  m <- match(c("x", "data", "subset", "weights"), names(mf), 0)
#  mf <- mf[c(1, m)]
#  mf$drop.unused.levels <- TRUE
#  mf[[1]] <- as.name("model.frame")
#  mf <- eval(mf, parent.frame())
#  y <- model.response(mf, "numeric")
#  w <- model.weights(mf)
#  mt <- attr(mf, "terms")
#  x <- x[[length(x)]]
#  if(length(x) == 3 && x[[1]] == as.name("|"))
#    x <- x[[2]]
#  x <- eval(substitute(x), mf)
#  curfit(x, y, w, ...)  
}

curfit.default <- function(x, y = NULL, w = NULL, s,
                           knots = NULL, n = 0,
                           from = min(x), to = max(x),
                           method = c("ls", "ss"), k = 3,
                           periodic = FALSE, ...) {
  method <- match.arg(method)
  if(method == "ls" && !missing(s)) {
    method <- "ss"
    warning("changing method to 'ss' since 's' was supplied")
  }
  iopt <- if(method == "ls") -1 else if(method == "ss") 0
  if(k - floor(k) > 0)
    warning("'k' must be an integer and will be truncated.")
  k <- as.integer(k)
  if(k < 1 || k > 5)
    stop("'k' must be between 1 and 5, inclusively")
  if(iopt == 0 && !is.null(knots)) iopt <- 1
  xlab <- deparse(substitute(x))
  ylab <- deparse(substitute(y))
  xyw <- xyw.coords(x, y, w)
  x <- xyw$x
  y <- xyw$y
  w <- xyw$w
  xin <- xyw$xin
  yin <- xyw$yin
  m <- length(x)
  dots <- list(...)
  if(from > min(x)) {
    warning("dropping values of x below ", from)
    keep <- x >= from
    x <- x[keep]
    y <- y[keep]
    w <- w[keep]
  }
  if(to < max(x)) {
    warning("dropping values of x above ", to)
    keep <- x <= to
    x <- x[keep]
    y <- y[keep]
    w <- w[keep]
  }
  if(missing(s)) {
    if(method == "ss") {
      stop("must supply 's' (smoothing parameter)")
    } else {
      s <- 0
    }
  }
  nest <- if("nest" %in% names(dots)) dots$nest else -1
  nest <- as.integer(max(nest, m + k + 1))
  if(iopt == -1) {
    if(missing(n) && is.null(knots))
      stop("must supply 'n' (number of knots) or 'knots' to use 'least squares splines'")
    xb <- min(x)
    xe <- max(x)
    if(!is.null(knots)) {
      if(!is.numeric(knots))
        stop("invalid argument for 'knots'")
      n <- length(knots)
      knots <- unique(sort(knots))
      eps <- .Machine$double.eps^0.5
      knots <- knots[knots > xb + eps & knots < xe - eps]
      knots <- if(periodic) {
        g <- length(knots)
        if(g >= k) {
          c(knots[rev(g + 1 - (1:k))] - xe + xb, xb, knots, xe, knots[1:k] + xe - xb)
        } else {
          stop("must have at least (", k, ") interior knots for periodic splines")
        }
      } else {
        c(rep(xb, k + 1), knots, rep(xe, k + 1))
      }
    } else {
      if(n < 0) stop("'n' must be greater than 0")
      interior <- seq(xb, xe, length = max(0, n - 2 * k))
      knots <- c(rep(xb, k), interior, rep(xe, k))
    }
    knots <- as.single(knots)
    n <- length(knots)
  } else if(iopt == 0) {
    if(!is.null(knots) || n > 0)
      warning("supplied 'knots' and/or 'n' ignored for smoothing splines")
    knots <- single(nest)
    n <- 0
  } else if(iopt == 1) {
    if(is.null(knots) || !is.numeric(knots) || length(knots) < 1)
      stop("invalid argument for 'knots'")
    n <- length(knots)
    knots <- as.single(knots)
  }
  coef <- single(nest)
  lwrk <- if("lwrk" %in% names(dots)) dots$lwrk else -1
  lwrk <- as.integer(max(lwrk, 2 * (m * (k + 1) + nest * (7 + 3 * k))))
  if(iopt == 1) {
    if(any(!c("wrk", "iwrk") %in% names(dots)))
      stop("Cannot update fit. Did you use 'update.dierckx'?")
    wrk <- as.single(c(dots$wrk[seq(n)], rep(0, lwrk - n)))
    iwrk <- as.integer(c(dots$iwrk[seq(n)], rep(0, nest - n)))
  } else {
    wrk <- single(lwrk)
    iwrk <- integer(nest)
  }
  val <- if(periodic) {
    ## percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
    .Fortran("percur",
             iopt = as.integer(iopt),
             m = as.integer(m),
             x = as.single(x),
             y = as.single(y),
             w = as.single(w),
             k = k,
             s = as.single(s),
             nest = as.integer(nest),
             n = as.integer(n),
             knots = knots,
             coef = coef,
             fp = single(1),
             wrk = wrk,
             lwrk = lwrk,
             iwrk = iwrk,
             ier = integer(1))
  } else {
    ## curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
    .Fortran("curfit",
             iopt = as.integer(iopt),
             m = as.integer(m),
             x = as.single(x),
             y = as.single(y),
             w = as.single(w),
             from = as.single(from),
             to = as.single(to),
             k = k,
             s = as.single(s),
             nest = as.integer(nest),
             n = as.integer(n),
             knots = knots,
             coef = coef,
             fp = single(1),
             wrk = wrk,
             lwrk = lwrk,
             iwrk = iwrk,
             ier = integer(1))
  }
  val$message <- if(val$ier == 0) {
    character(0)
  } else if(val$ier == -1) {
    "Spline returned is an interpolating spline."
  } else if(val$ier == -2) {
    msg1 <- sprintf("Spline returned is the weighted least-squares polynomial of degree %d.", k)
    msg2 <- sprintf("Upper bound on 's' is %f.", val$fp)
    paste(msg1, msg2)
  } else if(val$ier == 1) {
    stop("The required storage exceeds the available storage space specified by 'nest'.")
  } else if(val$ier == 2) {
    "A theoretically impossible result was found. possibly 's' is too small."
  } else if(val$ier == 3) {
    "The maximum number of iterations (20) has been  reached."
  } else if(val$ier == 10) {
    ""
  }
  if(val$ier == 10) {
    if(val$iopt > 1 || val$iopt < -1)
      val$message <- sprintf("%s Illegal option %d.", val$message, val$iopt)
    if(val$k < 1 || val$k > 5)
      val$message <- sprintf("%s Illegal 'k' (%d).", val$message, val$k)
    if(length(val$x) <= val$k)
      val$message <- sprintf("%s Length 'x' (%d) must be greater than 'k' (%d).",
                             val$message, m, val$k)
    if(val$nest < 2 * val$k + 2)
      val$message <- sprintf("%s Illegal value of 'nest' (%d).", val$message, val$nest)
    if(val$lwrk < with(val, (k + 1) * m + nest * (7 + 3 * k)))
      val$message <- sprintf("%s Illegal value of 'lwrk' (%d).", val$message, val$lwrk)
    if(!all(val$w > 0))
      val$message <- sprintf("%s Illegal weights.", val$message)
    if(val$iopt == -1) {
      n <- length(knots)
      if(n < 2 * k + 2 || n > min(val$nest, m + k + 1))
        val$message <- sprintf("%s Illegal number of knots (%d).", val$message, if(is.null(val$g)) 0 else val$g)
      if(val$knots[k + 2] < min(x) || val$knots[n - k - 1] > max(x) || all(x %in% val$knots))
        val$message <- sprintf("%s Illegal knot selection.", val$message)
    } else {
      if(s < 0)
        val$message <- sprintf("%s Illegal 's' (%f).", val$message, val$s)
      if(s == 0 && val$nest < m + val$k + 1)
        val$message <- sprintf("%s Illegal 'nest' (%d).", val$message, val$nest)
    }
    val$message <- sub("[ ]+", "", val$message)
  }
  if(val$ier > 0) stop(val$message)
  val$g <- with(val, length(knots[knots > min(x) & knots < max(x)]))
  sing <- names(unlist(sapply(val, attr, "Csingle")))
  val[sing] <- lapply(val[sing], as.double)
  val$method <- method
  val$periodic <- periodic
  val$routine <- "curfit.default"
  val$xlab <- xlab
  val$ylab <- ylab
  if(length(xin) > m) {
    val$x <- xin
    val$y <- yin
    val$sp <- predict.dierckx(val, newx = xin)
  } else {
    val$sp <- fitted(val)
  }
  class(val) <- "dierckx"
  val
}
