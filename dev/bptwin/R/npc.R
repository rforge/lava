###{{{ nonparcuminc

npc <- function(T,cause,same.cens=TRUE,sep=FALSE) {
  mtime <- apply(T[,1:2],1,max)
  ot <- order(mtime)
  mtime <- mtime[ot]
  T <- T[ot,]
  if (!sep) {
    time1 <- as.vector(T[,1:2]); status1 <- as.vector(T[,3:4])
    ud.cens1<-survfit(Surv(time1,status1==0)~+1);
    Gfit1<-cbind(ud.cens1$time,ud.cens1$surv)
    Gfit2 <- Gfit1<-rbind(c(0,1),Gfit1);
  } else {
    time1 <- as.vector(T[,1]); status1 <- as.vector(T[,3])
    ud.cens1<-survfit(Surv(time1,status1==0)~+1);
    time2 <- as.vector(T[,2]); status2 <- as.vector(T[,4])
    ud.cens2<-survfit(Surv(time2,status2==0)~+1);
    Gfit1<-cbind(ud.cens1$time,ud.cens1$surv)
    Gfit1<-rbind(c(0,1),Gfit1);
    Gfit2<-cbind(ud.cens2$time,ud.cens2$surv)
    Gfit2<-rbind(c(0,1),Gfit2);
  }
  cweights1<-fastapprox(Gfit1[,1],T[,1],Gfit1[,2])[[1]]
  cweights2<-fastapprox(Gfit2[,1],T[,2],Gfit2[,2])[[1]];
  weight11 <- apply(cbind(cweights1,cweights2),1,min)

  if (same.cens) {
    conc <- (T[,3]==cause[1])*(T[,4]==cause[2])/weight11
  } else {
    conc <-(T[,3]==cause[1])*(T[,4]==cause[2])/(cweights1*cweights2);
  }
  mtime <- mtime[!is.na(conc)]
  conc <- conc[!is.na(conc)]
  cbind(mtime,cumsum(conc)/length(conc))
}


nonparcuminc <- function(t,status,cens=0) {
  ord <- order(t); t <- t[ord]; status <- status[ord]
  ud.cens<-survfit(Surv(t,status==cens)~1)
  Gfit<-cbind(ud.cens$time,ud.cens$surv)
  Gfit<-rbind(c(0,1),Gfit);
  causes <- setdiff(unique(status),cens)
  cweight<-fastapprox(Gfit[,1],t,Gfit[,2])[[1]];
  cc <- t
  for (i in 1:length(causes)) {
    c1 <- status==causes[i]
    cc <- cbind(cc,cumsum(c1/cweight)/length(c1))
  }
  return(cc)
}

###}}} nonparcuminc
