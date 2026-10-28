#---include sources---
source("R/simulate.R")
source("R/kalman.R")


#---true values---
TT<-2000
phi<-0.8
sd_eta<- 1
sd_eps<- 0.5

#---simulator test---
set.seed(1)
res<-sim_ssm(TT,phi,sd_eta,sd_eps)
kf<-kalman_ll(res$y,phi,sd_eta,sd_eps)
print(sd(res$s)) #expected 1.6666
print(var(res$y)-var(res$s)) #expected 0.25
acf(res$s) #bars decaying 0.8^k

#----how is the filter doing? tracking----
plot(res$s,type="l")
lines(kf$a,col="red")
print(mean(abs(res$s-kf$a))) #0.37
print(mean(abs(res$s-res$y))) #0.415

#---innovations---
z<-(kf$v)/sqrt(kf$Fv)
print(mean(z)) #0
print(var(z)) #1
acf(z)



# ---variance path---
plot(kf$P,type="l")
print(tail(kf$P,1)) #0.204 expected

#---different dataset same variance path----
set.seed(3)
res2<-sim_ssm(TT,phi,sd_eta,sd_eps)
kf2<-kalman_ll(res2$y,phi,sd_eta,sd_eps)
print(all.equal(kf$P,kf2$P))

#---no randomness in kalman filter ---
kf3<-kalman_ll(res2$y,phi,sd_eta,sd_eps)
print(kf2$logL==kf3$logL)

