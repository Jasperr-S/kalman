#Efficient Importance Sampling
s<- seq(0,5,by=0.5)
w<- dnorm(3,s,0.5)*dnorm(s,0,1)/dnorm(s,2.4,sqrt(0.2))
print(w)
print(dnorm(3,0,sqrt(1.25)))