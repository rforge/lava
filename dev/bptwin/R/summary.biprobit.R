summary.biprobit <- function(object,level=0.05,...) {
  alpha <- level/2
  varcomp <- object$coef[length(coef(object)),1:2]
  varcomp <- rbind(object$model$tr(c(varcomp[1],varcomp[1]%x%cbind(1,1) + qnorm(1-alpha)*varcomp[2]%x%cbind(-1,1))))
  colnames(varcomp)[2:3] <- paste(c(alpha*100,100*(1-alpha)),"%",sep="")
  rownames(varcomp) <- ifelse(is.null(object$model$varcompname),"Variance component",object$model$varcompname)

  logit <- function(p) log(p/(1-p))
  tigol <- function(z) 1/(1+exp(-z))
  dlogit <- function(p) 1/(p*(1-p))
  probs <- function(p) {
    ##    S <- diag(2); S[1,2] <- S[2,1] <- exp(tail(p,1))
    S <- object$SigmaFun(p)
    m <- c(0,0)
    if (object$npar$intercept==1) m[1:2] <- p[1]
    if (object$npar$intercept==2) {
      m[1:2] <- p[c(1,1+object$npar$pred/2+1)]
    }
    mu.cond <- function(x) m[1]+S[1,2]/S[2,2]*(x-m[2])
    var.cond <- S[1,1]-S[1,2]^2/S[2,2]    
    conc <- pmvnorm(upper=m,sigma=S)
    marg <- pnorm(m[1],sd=S[1,1]^0.5)
    cond <- conc/marg
    logit(c(conc,cond,marg))
  }
  alpha <- level/2
  CIlab <- paste(c(alpha*100,100*(1-alpha)),"%",sep="")
  mycoef <- coef(object)
  prob <- probs(mycoef)
  Dprob <- jacobian(probs,mycoef)
  sprob <- diag((Dprob)%*%vcov(object)%*%t(Dprob))^0.5
  pp <- tigol(cbind(prob,prob-qnorm(1-alpha)*sprob,prob+qnorm(1-alpha)*sprob))
  rownames(pp) <- c("Concordance","Case-wise/Conditional","Marginal")
  colnames(pp) <- c("Estimate",CIlab)
  
  res <- list(object=object,varcomp=varcomp,prob=pp)
  class(res) <- "summary.biprobit"
  res
}
