###{{{ Dbvn

Dbvn <- function(p,design=function(p,...)
                 return(list(mu=cbind(p[1],p[1]),
                             dmu=cbind(1,1),
                             S=matrix(c(p[2],p[3],p[3],p[4]),ncol=2),
                             dS=rbind(c(1,0,0,0),c(0,1,1,0),c(0,0,0,1)))
                        ),                 
                        Y=cbind(0,0)) {
  mS <- design(p)
  U0 <- with(mS,.Call("biprobit0",
                          mu,
                          S,dS,Y,dmu,NULL,FALSE));
  return(c(U0,mS))
}

###}}} Dbvn
