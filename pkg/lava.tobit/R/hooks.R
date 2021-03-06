##' @export
lava.tobit.init.hook <- function(x,...) {
  x$attributes$binary <- list()
  return(x)
  ##  nodeDataDefaults(x,"binary") <- FALSE; x
}

##' @export
lava.tobit.sim.hook <- function(x,data,...) {  
  if (length(binary(x))>0)
    data[,binary(x)] <- (data[,binary(x)]>0)*1
  return(data)
}

##' @export
lava.tobit.color.hook <- function(x,subset=vars(x),...) {
  return(list(vars=intersect(subset,binary(x)),col="indianred1"))
}

##' @export
lava.tobit.estimate.hook <- function(x,data,weight,weight2,estimator,...) {
  dots <- list(...)
## Binary outcomes -> censored regression
  bin <- intersect(binary(x),vars(x))
  if (is.null(dim(data))) return(NULL)
  if (estimator%in%c("gaussian","tobit","normal")) {
    for (i in setdiff(endogenous(x),binary(x))) {
      if (is.character(data[,i]) | is.factor(data[,i])) { # Transform binary 'factor'
        y <- as.factor(data[,i])
        if (nlevels(y)==2) {
          data[,i] <- as.numeric(y)-1
          binary(x) <- i
        }
      }
    }
    if (length(bin)>0) {
      estimator <- "tobit"
      if (is.null(weight)) {        
        W <- data[,bin,drop=FALSE]; W[W==0] <- -1; colnames(W) <- bin
        weight <- lava.options()$threshold*W
      } else {
        ##        if (!all(binary(x)%in%colnames(data)))
        ##        W <- data[,binary(x),drop=FALSE]; W[W==0] <- -1; colnames(W) <- binary(x)
        ##        attributes(W)$weight2 <- weight
        ##        weight <- W
        ##          weight[,binary(x)] <- W
      }
      for (b in bin) {
        data[!is.na(data[,b]),b] <- 0
      }
      ##    data[,binary(x)] <- 0
      if (!is.null(weight2)) {
        estimator <- "tobitw"
      }
    }
  }
  
  ## Transform 'Surv' objects
  weight2 <- mynames <- NULL
  if (estimator%in%c("normal")) {
    for (i in setdiff(endogenous(x),bin)) {
      if (is.Surv(data[,i])) { 
        S <- data[,i]
        y1 <- S[,1]
        if (attributes(S)$type=="left")  {
          y2 <- y1
          y1[S[,2]==0] <- -Inf          
        }
        if (attributes(S)$type=="right") {
          y2 <- y1
          y2[S[,2]==0] <- Inf
        }
        if (attributes(S)$type=="interval2") {
          y2 <- S[,2]
        }
        if (!(attributes(S)$type%in%c("left","right","interval2"))) stop("Surv type not supported.")
        mynames <- c(mynames,i)
        weight2 <- cbind(weight2,y2)
        data[,i] <- y1
      }
    }
  }

  W <- NULL
  if (estimator%in%c("gaussian","tobit","tobitw")) {
    for (i in setdiff(endogenous(x),bin)) {
      if (is.Surv(data[,i])) { 
        estimator <- "tobit"
        S <- data[,i]
        y <- S[,1]
        if (attributes(S)$type=="left") 
          w <- S[,2]-1
        if (attributes(S)$type=="right") 
          w <- 1-S[,2]
        if (attributes(S)$type=="interval2") {
          w <- S[,3]; w[w==2] <- (-1)
        }
        mynames <- c(mynames,i)
        W <- cbind(W,w)
        data[,i] <- y
      }
    }
    if (length(W)>0) {
      colnames(W) <- mynames
      if (!is.null(weight)) {
        wW <- intersect(colnames(weight),colnames(W))
        if (length(wW)>0)
          weight[,wW] <- W[,wW]
        Wo <- setdiff(colnames(W),wW)
        if (length(Wo)>0)
        weight <- cbind(weight,W[,Wo,drop=FALSE])
      } else {
        weight <- W;
      }
    }
  }
  return(c(list(x=x,data=data,weight=weight,weight2=weight2,estimator=estimator),dots)) 
}
