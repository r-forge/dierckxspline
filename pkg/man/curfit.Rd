\name{curfit}
\alias{curfit}
\alias{percur}
\alias{curfit.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Curve fitting with splines}
\description{
General curve fitting with splines. Wrapper for function CURFIT.
}
\usage{
curfit(x, ...)
percur(x, ...)
\method{curfit}{default}(x, y = NULL, w = NULL, s,
       knots = NULL, n = 0, from = min(x), to = max(x),
       method = c("ls", "ss"), k = 3, periodic = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{data.frame}, \code{matrix}, or \code{numeric}
    vector. See details.}
  \item{y}{Optional numeric vector. See details.}
  \item{w}{Optional vector of weights}
  \item{s}{Smoothing parameter}
  \item{knots}{Optional vector of knots. See details.}
  \item{n}{The number of knots. Required only for method "ls".}
  \item{from}{Lower bound from which to fit spline.}
  \item{to}{Upper bound to which to fit spline.}
  \item{method}{Either "ss" for smoothing spline or "ls" for least
    squares spline. Partial matching is used to determine which method
    is supplied.}
  \item{k}{Degree of the spline. Valid options for \code{k} are 1 to 5,
    inclusively.}
  \item{periodic}{logical; if \code{TRUE} then \code{s(a) == s(b)}}
  \item{\dots}{Additional arguments used only in
    \code{update.curfit}. Otherwise, ignored.}
}
\details{
  As with \code{\link[stats]{smooth.spline}}, the \code{x} vector should
  contain at least four distinct values. \emph{Distinct} here means
  ``distinct after rounding to 6 significant digits'', i.e., \code{x}
  will be transformed to \code{unique(sort(signif(x, 6)))}, and \code{y}
  and \code{w} are pooled accordingly.
    
  For the default method, arguments \code{x} and \code{y} are supplied
  to \code{\link{xy.coords}} to determine abscissa and ordinate
  values. The actual function used is \code{DierckxSpline:::xyw.coords}
  which is not exported in the \code{NAMESPACE}.

  \code{curfit} determines a smooth spline approximation \code{s(x)}
  of degree \code{k}. See chapters 1 and 3 in the reference for definition of
  symbols.

  When supplying knots, the first \code{k+1} values should be
  \code{min(x)} and the last \code{n+k-1} value should be
  \code{max(x)}. See below for examples.

  \code{percur} is equivalent to calling \code{curfit(..., periodic=TRUE)}.
}
\value{
  An object of class \code{dierckx} with the following components:
  \item{iopt}{method used}
  \item{m}{length of 'x'}
  \item{x}{abscissa values}
  \item{y}{ordinate values}
  \item{w}{input weights}
  \item{from}{input \code{from} value}
  \item{to}{input \code{to} value}
  \item{s}{input smoothing parameter}
  \item{nest}{Estimated number of knots}
  \item{n}{Actual number of knots}
  \item{knots}{Knot locations.}
  \item{g}{Number of interior knots}
  \item{coef}{b-Spline coefficients. Use \code{coef.dierckx} to extract.}
  \item{fp}{sum of squares residuals. Use \code{deviance.dierckx} to
    extract.}
  \item{wrk}{Work space. Do NOT modify before call to
    \code{update.dierckx}}
  \item{lwrk}{Length of \code{wrk}. Do NOT modify before call to
    \code{update.dierckx}}
  \item{iwrk}{Integer work space. Do NOT modify before call to
    \code{update.dierckx}}
  \item{ier}{Error code. Should always be zero.}
  \item{method}{input \code{method} value}
  \item{routine}{Always 'curfit.default'}
  \item{xlab}{The x-label determined from \code{deparse(substitute(x)).}
  \item{ylab}{The y-label determined from \code{deparse(substitute(y)).}
}
\references{
Dierckx, P. (1991) \emph{Curve and Surface Fitting with Splines}, Oxford
Science Publications.
}
\author{Sundar Dorai-Raj}
\seealso{\code{\link{concon}}, \code{\link[stats]{spline}}, \code{\link[stats]{smooth.spline}}}
\examples{
x <- 0:24
y <- c(1.0,1.0,1.4,1.1,1.0,1.0,4.0,9.0,13.0,
       13.4,12.8,13.1,13.0,14.0,13.0,13.5,
       10.0,2.0,3.0,2.5,2.5,2.5,3.0,4.0,3.5)

z <- list()
for(k in c(3, 5, 2)) {
  z1 <- curfit(x, y, method = "ss", s = 1000, k = k)
  z2 <- update(z1, s = 60)
  z3 <- update(z2, s = 10)
  z4 <- update(z3, s = 30)
  z5 <- curfit(x, y, method = "ss", s = 30, k = k)
  z6 <- update(z5, s = 0)
  knots <- c(rep(0, k + 1), seq(3, 21, 3), rep(24, k + 1))
  z7 <- curfit(x, y, method = "ls", s = 30, knots = knots, k = k)
  z[[as.character(k)]] <- list(z1, z2, z3, z4, z5, z6, z7)
}

p <- unlist(z, recursive = FALSE)
n <- sapply(lapply(p, knots), length)
s <- sapply(p, "[[", "s")
i <- sapply(p, "[[", "iopt")
m <- ifelse(i == -1, "ls", ifelse(i == 0, "ss", "ss1"))
k <- sprintf("k = \%d", sapply(p, "[[", "k"))
g <- sprintf("\%s(s=\%d)", m, s, i)
sp <- data.frame(x = rep(x, times = length(p)),
                 y = rep(y, times = length(p)),
                 z = unlist(lapply(p, fitted)),
                 k = factor(rep(k, each = length(x))),
                 g = rep(g, each = length(x)))

library(lattice)
xyplot(z ~ x | k, data = sp, groups = g,
       panel = function(x, y, subscripts, groups, obs, ...) {
         panel.superpose(x, y, subscripts, groups, lwd = 3, type = "l", ...)
         x <- unique(x)
         y <- unique(obs)
         panel.xyplot(x, obs, pch = 16, cex = 1.2, col = "darkblue")
       },
       auto.key = list(space = "right", points = FALSE, lines = TRUE),
       obs = sp$y)

## periodic spline
set.seed(42)
n <- 100
r <- 1:n
x <- 0.01 * (r - 1)
e <- rnorm(n, 0, 0.1)
w <- rep(1/sd(e), n + 1)
y <- cos(2 * pi * x) + 0.25 * sin(8 * pi * x) + e
x <- c(x, 1)
y <- c(y, y[1])
kn <- seq(0.01, 0.99, length = 12)
f1 <- percur(x, y, w = w, method = "ss", s = 90, k = 5)

library(lattice)
top <- xyplot(y ~ x,
              panel = function(x, y, ...) {
                panel.abline(v = knots(f1), lty = 2, lwd = 3, col = "gray")
                panel.xyplot(x, y, pch = 16, col = "#800000", cex = 1.2)
                panel.xyplot(x, fitted(f1), type = "l", lwd = 3, col = "#000080")
              },
              par.settings = list(layout.widths = list(left.padding = 0, right.padding = 0)),
              scales = list(cex = 1.2),
              xlab = "", ylab = "")
newx <- seq(-2, 2, 0.01)
newy <- predict(f1, newx)
bot <- xyplot(newy ~ newx, type = "l",
              panel = function(...) {
                panel.abline(v = -2:2, lty = 2, col = "salmon", lwd = 3)
                panel.xyplot(...)
              },
              col = "#000080", lwd = 3,
              par.settings = list(layout.widths = list(left.padding = 0, right.padding = 0)),
              scales = list(cex = 1.2),
              xlab = "", ylab = "")
print(top, c(0, 0.2, 1, 1))
print(bot, c(0.008, 0, 0.992, 0.25), newpage = FALSE)

## example borrowed from ?smooth.spline
plot(cars$speed, cars$dist,
     main = "data(cars)  &  smoothing splines",
     xlab = "SPEED", ylab = "DISTANCE",
     cex.lab = 1.2, cex.axis = 1.2,
     cex.main = 2, cex = 1.5, col = "blue")
## This example has duplicate points, so avoid cv=TRUE
cars.spl.0 <- smooth.spline(cars$speed, cars$dist)
cars.spl.1 <- smooth.spline(cars$speed, cars$dist, df = 10)
cars.spl.2 <- curfit(cars$speed, cars$dist, s = 5e3)
newx <- seq(min(cars$speed), max(cars$speed), len = 200)
lines(predict(cars.spl.0, newx), col = "blue", lwd = 3, lty = 2)
lines(predict(cars.spl.1, newx), lty=2, col = "red", lwd = 3)
lines(newx, predict(cars.spl.2, newx), lty=2, col = "orange", lwd = 3)
legend(5, 120, c(paste("smooth.spline( * , df = ", round(cars.spl.1$df, 1), ")", sep = ""),
                 "smooth.spline( * , df = 10)", "curfit( * , s = 5e3)"),
       col = c("blue", "red", "orange"), lty = 2, lwd = 3,
       bg = 'bisque', cex = 1.5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
\keyword{optimize}
