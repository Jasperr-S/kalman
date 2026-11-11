source("R/simulate.R")
source("R/kalman.R")
source("R/particle.R")
set.seed(91)

TT<- 100
phi<- 0.8
sd_eta<- 0.3
sd_eps<- 0.5

#bench mark by kalmam
y<- sim_ssm(TT,phi,sd_eta,sd_eps)$y
k_ll<- kalman_ll(y,phi,sd_eta,sd_eps)$logL


#likelihood distribution and comparison (Jensen)
ptc_ll<- numeric(20)
for (i in 1:20){
  set.seed(i)
  ptc_ll[i]<- pf_ll(y,phi,sd_eta,sd_eps,10000)$logL
}

plot(1:20,ptc_ll,ylim=range(c(ptc_ll,k_ll)))
abline(h=k_ll)

#number of particles and variation
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

# -----The curves-----
# How much the phi hat moves if ll move by one sd?
# sd of ll about 0.1 for N=10000, about 0.74 for N=100

phi_grid<- seq(0.5,0.9,0.005)
kll_vec<- length(phi_grid)
for (i in seq_along(phi_grid)){
  kll_vec[i]<- kalman_ll(y,phi_grid[i],sd_eta,sd_eps)$logL
}
plot(phi_grid,kll_vec,type="l")
abline(h=max(kll_vec)-0.1,col="blue")
abline(h=max(kll_vec)-0.01,col="red")
gap_l<- range(phi_grid[kll_vec>max(kll_vec-0.01)])
gap_h<- range(phi_grid[kll_vec>max(kll_vec-0.1)])

#likelihood curve given by pf, fixed and varied seed
# varied curve should be half variance lower, about 0.29
pfll_f<- length(phi_grid)
pfll_v<- length(phi_grid)
for (i in seq_along(phi_grid)){
  set.seed(78)
  pfll_f[i]<-  pf_ll(y,phi_grid[i],sd_eta,sd_eps,100)$logL
}
for (i in seq_along(phi_grid)){
  pfll_v[i]<-  pf_ll(y,phi_grid[i],sd_eta,sd_eps,100)$logL
}
plot(phi_grid,kll_vec,type="l")
lines(phi_grid,pfll_f,col="blue")
lines(phi_grid,pfll_v,col="red")

# zoom in
phi_grid_fine<- seq(0.79,0.791,length.out=400)
pfll_z<- numeric(length(phi_grid_fine))
kll_z <- numeric(length(phi_grid_fine))
for (i in seq_along(phi_grid_fine)){
  set.seed(1)
  pfll_z[i]<-pf_ll(y[1:20],phi_grid_fine[i],sd_eta,sd_eps,100)$logL
  kll_z[i]<- kalman_ll(y[1:20],phi_grid_fine[i],sd_eta,sd_eps)$logL
}
plot(phi_grid_fine, pfll_z,type="l",col="blue",ylim=range(c(kll_z,pfll_z)))
lines(phi_grid_fine,kll_z,type="l")


#-----How the variation of estimation moves with N----

#---optimizer---
#kalman as the benchmark

neg_ll<- function(par,y){
  -kalman_ll(y,tanh(par),sd_eta,sd_eps)$logL
}

fit_ll<- optim(par=atanh(0.5),
               fn=neg_ll,
               y=y,method="BFGS",hessian = TRUE)

# variance of the estimation and the Hessain of ll
n<-500
Hll_vec<- numeric(n)
phi_vec<- numeric(n)
score_vec<- numeric(n)
TT<-100

for (i in 1:n){
  y_sim<- sim_ssm(TT,phi,sd_eta,sd_eps)$y
  fit_ll<- optim(par=atanh(0.5),
                 fn=neg_ll,
                 y=y_sim,method="BFGS",hessian=TRUE)
  phi_vec[i]<- tanh(fit_ll$par)
  Hll_vec[i]<- fit_ll$hessian[1,1]/(1-phi_vec[i]^2)^2
  score_vec[i] <- (kalman_ll(y_sim,phi+1e-4,sd_eta,sd_eps)$logL
      -kalman_ll(y_sim,phi-1e-4,sd_eta,sd_eps)$logL)/2e-4
}

print(mean(1/Hll_vec))
print(var(phi_vec))

print(var(score_vec))
print(mean(Hll_vec))



#optimizer for pf, with quasi newton mehod

N_grid<- c(100,300,1000,3000)

neg_llpf<- function(par,y,N,s){
  cat(tanh(par),"\n")
  set.seed(s)
  -pf_ll(y,tanh(par),sd_eta,sd_eps,N)$logL
}

failcount<-numeric(length(N_grid))
phi_mat<- matrix(NA,length(N_grid),50)
phi_sd<- numeric(length(N_grid))


for (i in seq_along(N_grid)){
  for (j in 1:50){
    fit_pf<- optim(par=atanh(0.5),
                   fn=neg_llpf,
                   y=y,method="BFGS",
                   N=N_grid[i],
                   s=j)
    failcount[i]<- failcount[i]+ (fit_pf$convergence!=0)
    phi_mat[i,j]<- fit_pf$par
  }

}
plot(tanh(phi_mat[4,1:50]))

phi_eff<- numeric(length(N_grid))
phi_var<- numeric(length(N_grid))
for (i in 1:length(N_grid)){
  phi_moved<- abs(tanh(phi_mat[i,])-0.5) >1e-2
  phi_eff[i]<- mean(phi_moved)
  phi_var[i]<- sd(tanh(phi_mat[i,phi_moved]))
}
plot(N_grid,phi_var,ylim=range(c(sd_eta/fit_ll$hessian ,phi_var)))
abline(h=sd_eta/fit_ll$hessian)



