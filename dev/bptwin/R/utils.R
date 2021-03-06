###{{{ Col
Col <- function (col, alpha = 0.2) 
{
    sapply(col, function(x) do.call(rgb, as.list(c(col2rgb(x)/255, 
        alpha))))
}
###}}} Col

###{{{ revdiag
revdiag <- function(x) {
  n <- ncol(x)
  x[cbind(rev(seq(n)),seq(n))]
}
"revdiag<-" <- function(x,...,value) {
  n <- ncol(x)
  x[cbind(rev(seq(n)),seq(n))] <- value
  x
}
###}}} revdiag

###{{{ ExMarg

ExMarg <- function(Y0,XX0,W0,dS0,midx1=seq(ncol(XX0)/2),midx2=seq(ncol(XX0)/2)+ncol(XX0)/2,eqmarg=TRUE,allmarg=FALSE) {  
  ii1 <- which(is.na(Y0[,2]) & !is.na(Y0[,1]))
  ii2 <- which(is.na(Y0[,1]) & !is.na(Y0[,2]))
  ii0 <- which(is.na(Y0[,1]) & is.na(Y0[,2]))
  margidx <- c(ii1,ii2)
  id1 <- id2 <-  NULL
  both <- setdiff(seq(nrow(Y0)),c(ii1,ii2,ii0))
  id <- seq_len(length(both))
  if (allmarg) {
    ##    id <- seq_len(length(both))
    id1 <- c(seq_len(length(ii1))+length(id), id)
    id2 <- c(seq_len(length(ii2))+length(id)+length(id1), id)
    ii1 <- c(ii1,both)
    ii2 <- c(ii2,both)
  }
  Y0_marg <- XX0_marg <- X0_marg1 <- X0_marg2 <- dS0_marg <- W0_marg <- NULL
  if (length(margidx)>0) {
    Y0_marg <- cbind(c(Y0[ii1,1],Y0[ii2,2]))
    X0_marg1 <- XX0[ii1,midx1,drop=FALSE]
    X0_marg2 <- XX0[ii2,midx2,drop=FALSE]
    dS0_marg <- dS0[,1,drop=FALSE]
    if (eqmarg) {
      XX0_marg <- rbind(X0_marg1,X0_marg2)
    } else {
      XX0_marg <- XX0[c(ii1,ii2),,drop=FALSE]
    }    
    if (!is.null(W0)) {
      W0_marg <- cbind(c(W0[ii1,1],W0[ii2,2]))
      W0 <- W0[-c(margidx,ii0),,drop=FALSE]
    }
    Y0 <- Y0[-c(margidx,ii0),,drop=FALSE]
    XX0 <- XX0[-c(margidx,ii0),,drop=FALSE]    
  }
  res <- list(Y0=Y0,XX0=XX0,W0=W0,
              Y0_marg=Y0_marg, XX0_marg=XX0_marg,
              X0_marg1=X0_marg1, X0_marg2=X0_marg2,
              dS0_marg=dS0_marg, W0_marg=W0_marg,
              id=id, idmarg=c(id1,id2),
              ii1=ii1,
              margidx=margidx)
}

###}}} ExMarg

###{{{ RoundMat

RoundMat <- function(cc,digits = max(3, getOption("digits") - 2),...) format(round(cc,max(1,digits)),digits=digits)

trMean <- function(b,blen) {
##  mytr <- function(x) x^2; dmytr <- function(x) 2*x
  mytr <- dmytr <- exp
  if (blen==0) return(b) 
  k <- length(b)
  Bidx <- seq_len(blen)+(k-blen)
  b[Bidx[1]] <- mytr(b[Bidx[1]])
  D <- diag(nrow=k)
  D[Bidx[1]:k,Bidx[1]] <- b[Bidx[1]]
  for (i in Bidx[-1]) {
    D[i:k,i] <- dmytr(b[i])
    b[i] <- b[i-1]+mytr(b[i])
  }
  attributes(b)$D <- D
  attributes(b)$idx <- Bidx
  return(b)
}

###}}} RoundMat

###{{{ multinomlogit

multinomlogit <- function(x) {
  n <- length(x)
  ex <- exp(x)
  sx <- sum(ex)+1
  f <- c(ex,1)
  df <- c(ex,0)
  res <- f/sx
  dg <- -ex/sx^2   
  gradient <- matrix(ncol=n,nrow=n+1)
  I <- diag(n+1)
  for (i in 1:(n)) {
    gradient[,i] <- df[i]*I[i,]/sx+dg[i]*f
  }
  attributes(res)$gradient <- gradient
  return(res)
}
###}}} multinomlogit

###{{{ Inverse

Inverse <- function(X,tol=1e-9) {
  n <- nrow(X)
  if (nrow(X)==1) {
    res <- 1/X
    if (det) attributes(res)$det <- X
    return(res)
  }
  svdX <- svd(X)
  id0 <- numeric(n)
  id0[svdX$d>tol] <- 1/svdX$d[svdX$d>tol]
  res <- with(svdX, v%*%diag(id0)%*%t(u))
  attributes(res)$sqrteig <- svdX$d
  attributes(res)$warning <- any(svdX$d<tol)
  return(res)
}

###}}} Inverse

###{{{ blockdiag
blockdiag <- function(x,...,pad=0) {
  if (is.list(x)) xx <- x  else xx <- list(x,...)
  rows <- unlist(lapply(xx,nrow))
  crows <- c(0,cumsum(rows))
  cols <- unlist(lapply(xx,ncol))
  ccols <- c(0,cumsum(cols))
  res <- matrix(pad,nrow=sum(rows),ncol=sum(cols))
  for (i in 1:length(xx)) {
    idx1 <- 1:rows[i]+crows[i]; idx2 <- 1:cols[i]+ccols[i]
    res[idx1,idx2] <- xx[[i]]
  }
  colnames(res) <- unlist(lapply(xx,colnames)); rownames(res) <- unlist(lapply(xx,rownames))
  return(res)
}
###}}} blockdiag

###{{{ Specials

Specials <- function(f,spec,split2="+",...) {
  tt <- terms(f,spec)
  pos <- attributes(tt)$specials[[spec]]
  if (is.null(pos)) return(NULL)
  x <- rownames(attributes(tt)$factors)[pos]
  st <- gsub(" ","",x)
  res <- unlist(strsplit(st,"[()]"))[2]
  if (is.null(split2)) return(res)
  unlist(strsplit(res,"+",fixed=TRUE))
}
###}}} Specials

###{{{ decomp.specials
decomp.specials <- function (x, pattern = "[()]", sep = ",", ...) 
  {
    st <- gsub(" ", "", x)
    if (!is.null(pattern)) 
      st <- rev(unlist(strsplit(st, pattern, ...)))[1]
    unlist(strsplit(st, sep, ...))
  }
###}}} decomp.specials
