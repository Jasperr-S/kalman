pf_ll<- function(y,phi,sd_eta,sd_eps,N){
  draws<- rnorm(N,0,sd_eta/sqrt(1-phi^2))
  loglik<- 0
  TT<-length(y)
  a_vec<- numeric(TT)
  P_vec<-numeric(TT)
  Fv_vec<- numeric(TT)
  v_vec<- numeric(TT)
  ESS_vec<-numeric(TT)
  var_ll_pred<- 0 #predicted variance of the loglik
  
  for (t in seq_along(y)){
    draws<- phi*draws+rnorm(N,0,sd_eta)
    a_pred<- mean(draws)
    P_pred<- var(draws)
    Fv_vec[t]<- P_pred+sd_eps^2
    
    v_vec[t]<- y[t]-a_pred
    
    w<- dnorm(y[t],draws,sd_eps)
    loglik<- loglik+log(sum(w)/N)
    w_bar<- w/sum(w)
    ESS_vec[t]<- 1/sum(w_bar^2)
    var_ll_pred<-var_ll_pred+1/ESS_vec[t]-1/N
    
    a_vec[t]<- sum(w_bar*draws)
    P_vec[t]<- sum(w_bar*(a_vec[t]-draws)^2)
    # draws<- sample(draws,N,replace=TRUE,prob=w_bar)
    #here use another resampler to maintain consistency of particles
    u<- runif(N)
    idx<- findInterval(u,cumsum(w_bar))+1
    idx<- pmin(idx,N)
    draws<- draws[idx]
    
  }
return(list(logL=loglik,a=a_vec,P=P_vec,v=v_vec
            ,Fv=Fv_vec,varlik=var_ll_pred,ESS=ESS_vec))
}
