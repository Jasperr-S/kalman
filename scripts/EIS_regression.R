lphi_fun<- function(s) dnorm(4,0,exp(s/2),log=TRUE)+dnorm(s,0,1,log=TRUE)
pyReal<- integrate(function(s) exp(lphi_fun(s)),-10,10)$value
N<- 10000

s<- seq(0,3,by=0.1)
fit<- lm(lphi_fun(s) ~ s+I(s^2))
plot(s,lphi_fun(s)); lines(s,fitted(fit))
sig2<- -1/2/coef(fit)[3]
mu<- coef(fit)[2]*sig2

#single draw
set.seed(1)
draws<- rnorm(N,mu,sqrt(sig2))
w<- exp(lphi_fun(draws)-dnorm(draws,mu,sqrt(sig2),log=TRUE))
pEST<- mean(w)
ESS<- sum(w)^2/sum(w^2)

#bootstrap as benchmark
set.seed(1)
draws2<- rnorm(N,0,1)
w2<- dnorm(4,draws2,exp(draws2/2))
ESS2<- sum(w2)^2/sum(w2^2)

#let g choose its own points, with common random numbers
set.seed(3)
u<-rnorm(100,0,1)
muG<- 0
sig2G<- 1
TT<- 10
muVec<- numeric(TT)
sig2Vec<- numeric(TT)
for (i in 1:TT){
  drawsG<- u*sqrt(sig2G)+muG
  #drawsG<-rnorm(100,muG,sqrt(sig2G))
  fitG<- lm(lphi_fun(drawsG)~drawsG+I(drawsG^2))
  sig2G<- -1/2/coef(fitG)[3]
  muG<- sig2G*coef(fitG)[2]
  muVec[i]<- muG
  sig2Vec[i]<- sig2G
  if (sig2G<0) stop("parabola opens upward")
}
draws3<- rnorm(N,muG,sqrt(sig2G))
w3<- exp(lphi_fun(draws3)-dnorm(draws3,muG,sqrt(sig2G),log=TRUE))
pEST3<- mean(w3)
ESS3<- sum(w3)^2/sum(w3^2)


#smoothness check
set.seed(78)
u<-rnorm(100,0,1)
u2<- rnorm(1000,0,1)
lphi_funE<- function(s,beta) dnorm(4,0,beta*exp(s/2),log=T)+dnorm(s,0,1,log=T)
eis_one<- function(beta,u,u2){
  sig2E<- 1
  muE<- 0
  for (i in 1:10){
    drawsE<- u*sqrt(sig2E)+muE
    fitE<-lm(lphi_funE(drawsE,beta)~drawsE+I(drawsE^2))
    sig2E<- -1/2/coef(fitE)[3]
    muE<- sig2E*coef(fitE)[2]
    if(sig2E<0) stop ("parabola opens upward")
  }
  drawsE2<- u2*sqrt(sig2E)+muE
  w<- exp(lphi_funE(drawsE2,beta)-dnorm(drawsE2,muE,sqrt(sig2E),log=T))
  return(log(mean(w)))
}

betas<- seq(0.5,1,by=0.01)
ll_CRN<-numeric(length(betas)) #common random numbers
ll_FD<-numeric(length(betas)) #fresh draws
ll_B<- numeric(length(betas)) #bootstrap with fixed draws
ll_T<- numeric(length(betas))
for (b in seq_along(betas)){
  beta<- betas[b]
  ll_CRN[b]<- eis_one(beta,u,u2)
  ll_T[b]<- log(integrate(function(s) exp(lphi_funE(s,beta=beta)),-10,10)$value)
  ll_B[b]<- log(mean(dnorm(4,0,beta*exp(u2/2))))
}

for (b in seq_along(betas)){
  beta<- betas[b]
  u<-rnorm(100,0,1)
  u2<-rnorm(1000,0,1)
  ll_FD[b]<- eis_one(beta,u,u2)
}

plot(betas,ll_CRN-ll_T,type="l",ylim=range(c(ll_CRN-ll_T,ll_B-ll_T,ll_FD-ll_T)),ylab="gap")
lines(betas,ll_B-ll_T,type="l",col="red")
lines(betas,ll_FD-ll_T,type="l",col="blue")
abline(h=0,col="green")







