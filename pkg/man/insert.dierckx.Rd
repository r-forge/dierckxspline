\name{insert.dierckx}
\alias{insert}
\alias{insert.dierckx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spline Knot Insertion}
\description{
  Inserts a knot into a spline object.
}
\usage{
insert(object, ...)
\method{insert}{dierckx}(object, at, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{dierckx}.}
  \item{at}{A vector of knots to insert.}
  \item{\dots}{ignored}
}
\details{
}
\value{
  An updated spline representation.
}
\references{
Dierckx, P. (1991) \emph{Curve and Surface Fitting with Splines}, Oxford
Science Publications.
}
\author{Sundar Dorai-Raj}
\seealso{\code{\link{curfit}}}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
\keyword{math}
