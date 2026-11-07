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

plot(1:20,ptc_ll)
abline(v=k_ll)