pf_ll<- function(y,phi,sd_eta,sd_eps,N){
  draws<- rnorm(N,0,sd_eta/sqrt(1-phi^2))
  loglik<- 0
  TT<-length(y)
  a_vec<- numeric(TT)
  P_vec<-numeric(TT)
  Fv_vec<- numeric(TT)
  v_vec<- numeric(TT)
  
  for (t in seq_along(y)){
    draws<- phi*draws+rnorm(N,0,sd_eta)
    a_pred<- mean(draws)
    P_pred<- var(draws)
    Fv_vec[t]<- P_pred+sd_eps^2
    
    v_vec[t]<- y[t]-a_pred
    
    w<- dnorm(y[t],draws,sd_eps)
    loglik<- loglik+log(sum(w)/N)
    w_bar<- w/sum(w)
    a_vec[t]<- sum(w_bar*draws)
    P_vec[t]<- sum(w_bar*(a_vec[t]-draws)^2)
    draws<- sample(draws,N,replace=TRUE,prob=w_bar)
    
  }
return(list(logL=loglik,a=a_vec,P=P_vec,v=v_vec,Fv=Fv_vec))
}
