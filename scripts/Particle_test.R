draws<-rnorm(10000,2,0.5)
print(mean(draws))
print(var(draws))
print(mean(draws>3))
hist(draws,breaks=50)

#---predict---
phi<-0.9
sig_eta<-0.3
new<-phi*draws+rnorm(10000,0,sig_eta)
print(mean(new))
print(var(new))
hist(new,breaks=50)

new2<-sin(3*draws)+rnorm(10000,0,sig_eta)
hist(new2,breaks=50)

#---weights and update, y_t+1=2.5---
y<-4
sig_eps<-0.5
w<-dnorm(y,new,sig_eps)
wbar<-w/sum(w)
m_upd<-sum(wbar*new)
var_upd<-sum(wbar*(new-m_upd)^2)
print(c(m_upd,var_upd))

#---what if kalman?---
P<-phi^2*0.25+sig_eta^2
R<-sig_eps^2
K<-P/(P+R)
m_K <-1.8+K*(y-1.8)
var_K<- P*(1-K)
print(c(m_K,var_K))


#---ESS effective sample size----
y_list<-seq(0,4,0.05)
ESS_list<- numeric(length(y_list))
ESS_analy<-numeric(length(y_list))
sd_mean<-numeric(length(y_list))
for (i in seq_along(y_list)){
  w<- dnorm(y_list[i],new,sig_eps)
  wbar<-w/sum(w)
  ESS_list[i]<-1/sum(wbar^2)
  m_upd<-sum(wbar*new)
  var_upd<-sum(wbar*(new-m_upd)^2)
  sd_mean[i]<-sqrt(var_upd/ESS_list[i])
  ESS_analy[i]<-10000*R/(P+R)*sqrt((2*P+R)/R)*exp(-P*(y_list[i]-mean(new))^2/(P+R)/(2*P+R))
}
plot(y_list,sd_mean,type="l")
plot(y_list,ESS_list,type="l",ylim = range(c(ESS_list,ESS_analy)))
lines(y_list,ESS_analy,type="l",col="red")

  
  
  
  

      