\name{mixture}
\alias{mixture}
\title{Estimate mixture latent variable model}
\description{
  Estimate parameters in a mixture of latent variable models via the EM algorithm.
}
\usage{
mixture(x, data, k=length(x), control, FUN, type=c("standard","CEM","SEM"),...)
}
\arguments{
  \item{x}{List of \code{lvm} objects. If only a single \code{lvm}
    object is given, then a \code{k}-mixture of this model is fitted
    (free parameters varying between mixture components)}
  \item{data}{\code{data.frame}}
  \item{k}{Number of mixture components}
  \item{control}{Optimization parameters (see details)}
  \item{FUN}{See details below}
  \item{type}{Type of EM algorithm (standard, classification, stochastic)}
  \item{\dots}{Additional arguments parsed to lower-level functions}
}
\value{
  A \code{mixture} object
}
\details{
The performance of the EM algorithm can be tuned via the \code{control}
argument, a list where a subset of the following members can be altered:

\describe{
\item{start}{Optional starting values}
\item{nstart}{Evaluate \code{nstart} different starting values and run the
  EM-algorithm on the parameters with largest likelihood}
\item{tol}{Convergence tolerance of the EM-algorithm.  The algorithm is
  stopped when the absolute change in likelihood and parameter (2-norm)
  between successive iterations is less than \code{tol}}
\item{iter.max}{Maximum number of iterations of the EM-algorithm}
\item{gamma}{Scale-down (i.e. number between 0 and 1) of the step-size
  of the Newton-Raphson algorithm in the M-step}
\item{trace}{Trace information on the EM-algorithm is printed on every
  \code{trace}th iteration}
}

Note that the algorithm can be aborted any time (C-c) and still be saved
(via on.exit call).
}
\examples{
set.seed(1)
m0 <- lvm(list(y~x+z,x~z))
distribution(m0,~z) <- binomial.lvm()
d <- sim(m0,500,p=c("y<-z"=2,"y<-x"=1))

## unmeasured confounder example
m <- baptize(lvm(y~x));
covariance(m,~x) <- "v"
intercept(m,~x+y) <- NA

\donttest{
M <- mixture(m,k=2,data=d,control=list(trace=1,tol=1e-4))
summary(M)
lm(y~x,d)
## True slope := 1
}
}
\seealso{\code{mvnmix}}
\author{Klaus K. Holst}
\keyword{models}
\keyword{regression}
