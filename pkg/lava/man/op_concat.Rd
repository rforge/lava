% Generated by roxygen2 (4.0.1): do not edit by hand
\name{\%+\%}
\alias{\%++\%}
\alias{\%+\%}
\title{Concatenation operator}
\usage{
x \%+\% y
}
\arguments{
\item{x}{First object}

\item{y}{Second object of same class}
}
\description{
For matrices a block-diagonal matrix is created. For all other
data types he operator is a wrapper of \code{paste}.
}
\details{
Concatenation operator
}
\examples{
matrix(rnorm(25),5)\%+\%matrix(rnorm(25),5)
"Hello "\%+\%" World"
}
\author{
Klaus K. Holst
}
\seealso{
\code{blockdiag}, \code{\link{paste}}, \code{\link{cat}},
}
\keyword{misc}
\keyword{utilities}

