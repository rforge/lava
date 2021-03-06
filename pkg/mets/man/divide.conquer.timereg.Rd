% Generated by roxygen2 (4.0.1): do not edit by hand
\name{divide.conquer.timereg}
\alias{divide.conquer.timereg}
\title{Split a data set and run function from timereg and aggregate}
\usage{
divide.conquer.timereg(func = NULL, data, size, ...)
}
\arguments{
\item{func}{called function}

\item{data}{data-frame}

\item{size}{size of splits}

\item{...}{Additional arguments to lower level functions}
}
\description{
Split a data set and run function of cox-aalen type and aggregate results
}
\examples{
library(timereg)
data(TRACE)
\donttest{
a <- divide.conquer.timereg(prop.odds,TRACE,
	     formula=Surv(time,status==9)~chf+vf+age,n.sim=0,size=200)
coef(a)
a2 <- divide.conquer.timereg(prop.odds,TRACE,
	     formula=Surv(time,status==9)~chf+vf+age,n.sim=0,size=500)
coef(a2)
par(mfrow=c(1,1))
plot(a,xlim=c(0,8),ylim=c(0,0.01))
par(new=TRUE)
plot(a2,xlim=c(0,8),ylim=c(0,0.01))
}
}
\author{
Thomas Scheike, Klaus K. Holst
}

