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

# -----The curves
# How much the real value moves if ll move by one sd?
# sd=0.1 for N=10000, about 0.74 for N=100
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







