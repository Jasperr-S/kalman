source("R/simulate.R")
source("R/kalman.R")
dir.create("figures", showWarnings = FALSE)

#---true parameters---
phi_true<-0.8
TT<-2000
sd_eta_true<- 1
sd_eps_true<- 0.5
res<- sim_ssm(TT,phi_true,sd_eta_true,sd_eps_true)

#---verify wrt each parameter---
phi_grid<-seq(from=0.5,to=0.95,by=0.005)
ll<-numeric(length(phi_grid))
pen<-ll
fitq<-ll

for (i in seq_along(phi_grid)){
  k<-kalman_ll(res$y,phi_grid[i],sd_eta_true,sd_eps_true)
  ll[i]<- k$logL
  pen[i]<- -sum(log(k$Fv))/2
  fitq[i]<- -0.5*sum(k$v^2/k$Fv)
}

png("figures/loglik-phi.png",width=800,height=600)
plot(phi_grid,ll,type="l",xlab="phi",ylab="log-likelihood")
abline(v=phi_true,col="red",lty=2)
abline(v=phi_grid[which.max(ll)],col="blue")
dev.off()

#----log-likelihood decomposition----

plot(phi_grid,y=fitq,type="l",col="blue",ylim=range(c(fitq,pen))
     ,xlab="phi",ylab="contribution to logL")
lines(phi_grid,pen,col="red")
legend("bottomleft", legend = c("fit term", "penalty term"),
       col = c("blue", "red"), lty = 1)

#----first order difference decomposition----
d_ll<-diff(ll)/diff(phi_grid)
d_pen<-diff(pen)/diff(phi_grid)
d_fitq<-diff(fitq)/diff(phi_grid)
phi_mid <- phi_grid[-1]-diff(phi_grid)/2

png("figures/loglik_FOC_decomp.png",width=800,height = 600)
plot(phi_mid,d_fitq,type="l",col="blue",
     ylim=range(c(d_fitq,d_ll,d_pen)),
     xlab="phi",ylab="d(logL)/d(phi)")
lines(phi_mid,d_pen,col="red")
lines(phi_mid,d_ll,col="black")
abline(h=0,lty=3)
abline(v=phi_true,col="grey",lty=2)
legend("bottomleft", legend = c("fit term", "penalty term"),
       col = c("blue", "red"), lty = 1)
dev.off()

#----convergence rate visualization----
TT_small<-200
TT_large<-20000
ll_large<-numeric(length(phi_grid))
ll_small<-numeric(length(phi_grid))
res_large<- sim_ssm(TT_large,phi_true,sd_eta_true,sd_eps_true)
res_small<- sim_ssm(TT_small,phi_true,sd_eta_true,sd_eps_true)
for (i in seq_along(phi_grid)){
  k_large<-kalman_ll(res_large$y,phi_grid[i],sd_eta_true,sd_eps_true)
  ll_large[i]<- k_large$logL
  k_small<-kalman_ll(res_small$y,phi_grid[i],sd_eta_true,sd_eps_true)
  ll_small[i]<- k_small$logL
}

plot(phi_grid,ll,type="l",xlab="phi",ylab="log-likelihood"
     ,ylim=range(c(ll,ll_large,ll_small)))
lines(phi_grid,ll_large,col="blue")
lines(phi_grid,ll_small,col="red")



png("figures/convergence_rate.png",width=800,height=600)
plot(phi_grid, ll_small - max(ll_small), type = "l", col = "red",
     ylim = c(-50, 2), xlab = "phi", ylab = "logL (centred at max)")
lines(phi_grid, ll       - max(ll),       col = "black")
lines(phi_grid, ll_large - max(ll_large), col = "blue")
abline(v = phi_true, lty = 2, col = "grey")
legend("bottomright", legend = c("T=200", "T=2000", "T=20000"),
       col = c("red", "black", "blue"), lty = 1)
abline(h=-1.92,lty=3)
dev.off()

#----sd-eta----
sd_eta_grid<- seq(from=0.5,to=2,by=0.005)
ll_eta<- numeric(length(sd_eta_grid))
for (i in seq_along(sd_eta_grid)){
  ll_eta[i]<-kalman_ll(res$y,phi_true,sd_eta_grid[i],sd_eps_true)$logL
}
png("figures/loglik-sd_eta.png",width=800,height=600)
plot(sd_eta_grid,ll_eta,type="l",xlab="sd_eta",ylab="log_likelyhood")
abline(v=sd_eta_true,col="red",lty=2)
abline(v=sd_eta_grid[which.max(ll_eta)],col="blue")
dev.off()


#---optimizer---

neg_ll<- function(par,y){
  -kalman_ll(y,tanh(par[1]),exp(par[2]),exp(par[3]))$logL
  
}

fit<- optim(par=c(atanh(0.5),log(2),log(2)),
            fn=neg_ll,
            y=res$y)

print(c(tanh(fit$par[1]),exp(fit$par[2]),exp(fit$par[3])))
print(fit$convergence)






