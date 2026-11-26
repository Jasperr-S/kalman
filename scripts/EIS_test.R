#Efficient Importance Sampling

#the weight
s<- seq(0,5,by=0.5)
w<- dnorm(3,s,0.5)*dnorm(s,0,1)/dnorm(s,2.4,sqrt(0.2))
print(w)
print(dnorm(3,0,sqrt(1.25)))

# the regression for a bell curve
sFine<- seq(0,5,by=0.1)
#log(p(y|s)*p(s|y_t-1))
lphi<- dnorm(3,sFine,0.5,log=TRUE)+dnorm(sFine,0,1,log=TRUE)
fit<- lm(lphi ~ sFine+I(sFine^2))
sig2<- -1/(2*coef(fit)[3])
mu<- sig2*coef(fit)[2]
cat("predicted sigma^2=",sig2)
cat("predicted miu=",-coef(fit)[2]/(2*coef(fit)[3]))
summary(fit)$r.squared
cat("p(yt|yt-1)=",exp(coef(fit)[1] + 0.5*log(2*pi*sig2) + mu^2/(2*sig2)))

