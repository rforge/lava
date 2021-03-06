ClaytonOakes <- function(formula,data=parent.frame(),id,var.formula=~1,cuts=NULL,type="co",start,control=list(),...) {
  
  mycall <- match.call()
  formulaId <- Specials(formula,"id")
  
  formulaStrata <- Specials(formula,"strata")
  formulaSt <- paste("~.-id(",formulaId,")",
                     "-strata(",paste(formulaStrata,collapse="+"),")")
  formula <- update(formula,formulaSt)
  if (!is.null(formulaId)) {
    decompId <- decomp.specials(formulaId)
    if (length(decompId)>1) {
      var.formula <- as.formula(decompId[1])
      formulaId <- decompId[2]
    }    
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) {
    strata <- formulaStrata
    mycall$strata <- strata
  }
  if (missing(id)) stop("Missing 'id' variable")
  
  timevar <- terms(formula)[[2]]
  if (is.call(timevar)) {
    delayedentry <- (length(timevar)==4)*1
    entry <- NULL
    if (delayedentry==1)
      entry <- as.character(timevar[[2]])
    causes <- timevar[[3+delayedentry]]
    timevar <- timevar[[2+delayedentry]]
  }  
  timevar <- as.character(timevar)
  causes <- as.character(causes)
  covars <- as.character(attributes(terms(formula))$variables)[-(1:2)]
  X <- NULL
  nbeta <- 0
  if (length(covars)>0) {
    X <- model.matrix(as.formula(paste("~-1+",covars)),data)
    nbeta <- ncol(X)
  }
  ngamma <- 0
  Z <- model.matrix(var.formula,data)
  ngamma <- ncol(Z)
  
  if (is.data.frame(data)) {
    mydata <- data.frame(T=data[,timevar],status=data[,causes],cluster=data[,id],entry=0)
    if (!is.null(entry)) {
      mydata$entry <- data[,entry]
    }
  } else {
    mydata <- data.frame(T=get(timevar,env=data),status=get(causes,env=data),cluster=get(id,env=data),entry=0)
    if (!is.null(entry))      
      mydata$entry <- get(entry,env=data)
  }
  if (is.null(cuts)) {
    cuts <- c(0,max(mydata$T))
  }
  if (max(mydata$T)>tail(cuts,1)) stop("Interval does not embed time observations")
  if (any(with(mydata, T<entry))) stop("Entry time occuring after event")

  ucluster <- unique(mydata$cluster)
  
  npar <- length(cuts)-1
  if (!is.null(X)) npar <- npar+ncol(X)
  npar <- npar+ncol(Z)
  p0 <- rep(0,npar)
  if (!missing(start)) p0 <- c(start,rep(0,max(0,length(npar)-length(start))))
  
  obj <- function(p) {
    varpar <- p[seq(ngamma)]
    p <- p[-seq(ngamma)]
    ##    theta0 <- rep(exp(varpar),length(ucluster));
    theta0 <- exp(Z%*%varpar)
    multhaz <- rep(1,nrow(mydata))
    if (!is.null(X)) {
      nbeta <- ncol(X)
      beta <- p[seq(nbeta)]
      p <- p[-seq(nbeta)]
      multhaz <- exp(X%*%beta)
    }
    res <- .Call("claytonoakes",
           ds=mydata$status,ts=mydata$T,es=mydata$entry,
           allcs=mydata$cluster,cs=ucluster,
           cuts=cuts,hs=exp(p),mulths=multhaz,
           var=theta0,DUP=FALSE)$logLik
    return(-res)
  }

  tryCatch(opt <- nlminb(p0,obj,control=control),error=function(x) NA)
  I <- hessian(obj,opt$par)
  ee <- eigen(I);
  threshold <- 1e-12
  idx <- ee$values>threshold
  ee$values[idx] <- 1/ee$values[idx];
  if (!all(idx))
    ee$values[!idx] <- 0
  V <- with(ee, vectors%*%diag(values)%*%t(vectors))
  res <- list(coef=opt$par,vcov=V,cuts=cuts,nbeta=nbeta,ngamma=ngamma,betanames=colnames(X),gammanames=colnames(Z),opt=opt)
  class(res) <- "claytonoakes"
  return(res)
}

##################################################

print.claytonoakes <- function(x,...) {
  print(summary(x))
}

print.summary.claytonoakes <- function(x,...) {
  print(x$coef[,c(1,3,4)])
}

summary.claytonoakes <- function(obj,...) {
  mycoef <- matrix(nrow=length(obj$coef),ncol=4)
  mycoef[,1:2] <- cbind(obj$coef,sqrt(diag(obj$vcov)))
  mycoef[,3:4] <- cbind(mycoef[,1]-qnorm(0.975)*mycoef[,2],mycoef[,1]+qnorm(0.975)*mycoef[,2])
  colnames(mycoef) <- c("Estimate","Std.Err","2.5%","97.5%")
  if (length(obj$cuts))
  cutnames <- levels(cut(0,breaks=obj$cuts))
  rownames(mycoef) <- c(paste("log-Var:",obj$gammanames,sep=""),obj$betanames,cutnames)
  mycoef[-seq(obj$ngamma),] <- exp(mycoef[-seq(obj$ngamma),])
  res <- list(coef=mycoef)
  class(res) <- "summary.claytonoakes"
  res
}

plot.claytonoakes <- function(x,chaz=TRUE,add=!is.null(dev.list()),col="darkblue",...) {
  haz <- summary(x)$coef[-seq(x$nbeta+x$ngamma),,drop=FALSE]
  t <- x$cuts
  L <- approxfun(t,f=1,cumsum(c(0,haz[,1]*diff(t))),method="linear")
  if (add) {
    lines(t,L(t),col=col,...)
  } else {
    plot(t,L(t),type="l",col=col,...)
  }
  invisible(x)  
}
