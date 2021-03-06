% Generated by roxygen2 (4.0.1): do not edit by hand
\name{lava.options}
\alias{lava.options}
\title{Set global options for \code{lava}}
\usage{
lava.options(...)
}
\arguments{
\item{\dots}{Arguments}
}
\value{
\code{list} of parameters
}
\description{
Extract and set global parameters of \code{lava}. In particular optimization
parameters for the \code{estimate} function.
}
\details{
\itemize{ \item \code{param}: 'relative' (factor loading and variance of one
endogenous variables in each measurement model are fixed to one), 'absolute'
(mean and variance of latent variables are set to 0 and 1, respectively),
'hybrid' (intercept of latent variables is fixed to 0, and factor loading of
at least one endogenous variable in each measurement model is fixed to 1),
'none' (no constraints are added) \item \code{silent}: Set to \code{FALSE}
to disable various output messages \item ...  } see \code{control} parameter
of the \code{estimate} function.
}
\examples{
\dontrun{
lava.options(iter.max=100,silent=TRUE)
}
}
\author{
Klaus K. Holst
}
\keyword{models}

