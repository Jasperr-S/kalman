# simulate state space model

sim_ssm <- function (TT,phi,sigma_eta,sigma_eps){
  s<-numeric(TT)
  y<-numeric(TT)
  s0<-rnorm(1,0,sigma_eta/sqrt(1-phi^2))
  s[1]<-s0*phi+rnorm(1,0,sigma_eta)
  y[1]<-s[1]+rnorm(1,0,sigma_eps)
  for (t in 2:TT){
    s[t]<-s[t-1]*phi+rnorm(1,0,sigma_eta)
    y[t]<-s[t]+rnorm(1,0,sigma_eps)
    
  }
  return(list(y=y,s=s))
}