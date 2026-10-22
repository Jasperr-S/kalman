kalman_ll<-function(y,phi,sd_eta,sd_eps){
  TT<-length(y)
  logL<- -log(2*pi)*TT/2
  a_vec<-numeric(TT)
  Fv_vec<-numeric(TT)
  P_vec<-numeric(TT)
  v_vec<-numeric(TT)
  
  
  a<-0 #initialize a_(1|0)
  P<-phi^2*sd_eta^2/(1-phi^2)+sd_eta^2
  
  
  for (i in 1:TT){
    v_vec[i]<-y[i]-a
    Fv_vec[i]<-P+sd_eps^2
    K<-P/Fv_vec[i]
    a_vec[i]<-a+K*v_vec[i]
    P_vec[i]<-P*(1-K)
    logL<-logL-log(Fv_vec[i])/2-v_vec[i]^2/Fv_vec[i]/2
    a<-phi*a_vec[i]
    P<-phi^2*P_vec[i]+sd_eta^2
  }
  
  return(list(logL=logL,a=a_vec,P=P_vec,v=v_vec,Fv=Fv_vec))
  
  
  
}