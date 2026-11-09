source("R/simulate.R")
source("R/kalman.R")
source("R/particle.R")
set.seed(91)

TT<- 100
phi<- 0.8
sd_eta<- 0.3
sd_eps<- 0.5

y<- sim_ssm(TT,phi,sd_eta,sd_eps)$y
k_ll<- kalman_ll(y,phi,sd_eta,sd_eps)$logL

ptc_ll<- numeric(20)
for (i in 1:20){
  set.seed(i)
  ptc_ll[i]<- pf_ll(y,phi,sd_eta,sd_eps,10000)$logL
}

plot(1:20,ptc_ll,ylim=range(c(ptc_ll,k_ll)))
abline(h=k_ll)


ptc_ll_ss<- numeric(20)
pred_ll<-numeric(20)
for (i in 1:20){
  set.seed(i)
  out<- pf_ll(y,phi,sd_eta,sd_eps,100)
  ptc_ll_ss[i]< -out$logL
  pred_ll[i]<- out$varlik
}
print(mean(pred_ll))
print(var(ptc_ll_ss))


print(sd(ptc_ll_ss)/sd(ptc_ll))
