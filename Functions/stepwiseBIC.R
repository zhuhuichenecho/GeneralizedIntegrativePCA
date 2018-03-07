GIPCABIC<-function(data,rankj,ranka,tol){
  
  
  
  N = nrow(do.call(cbind,data)); P = ncol(do.call(cbind,data)); K = rankj; 
  
  for(i in 1:length(data)){
    
    eval(parse(text = paste("K.ind",i,"=ranka[",i,"]",sep = "")))
  }
  
  D = sapply(data,ncol)
  dimention = c(0,cumsum(D))
  
  for (d in 1:(length(dimention)-1)){
    
    eval(parse(text = paste("N.ind",d,"=nrow(na.omit(X[,",(dimention[d]+1),":",dimention[d+1],"]))",sep = "")))
    
  }
  
  result = EPCAJIVEMissbio(data,rankj,ranka,D = D,family = family, tol = tol,max.iter =250,n.size = n.size)
  
  
  final.BIC.EPCAJIVE2 = -2*result$loglik+log(N.ind1*D[1]+N.ind2*D[2])*
    (P+ifelse(K==0,0,sum((N-K):(N-1))+P*(K))+
       ifelse(K.ind1==0,0,sum((N.ind1-(K.ind1+K)):(N.ind1-K-1))+(K.ind1)*P/2)+
       ifelse(K.ind2==0,0,sum((N.ind2-(K.ind2+K)):(N.ind2-K-1))+(K.ind2)*P/2))
  
  
  return(final.BIC.EPCAJIVE2)
  
}
