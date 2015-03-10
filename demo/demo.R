

##simulation study
#training set
a1 = 0.9; a2 = -0.8; f = 0.1
dat = c()
dat[1] = rnorm(1);dat[2] = rnorm(1)
for (k in 2:1000) {
  dat[k+1] = a1*dat[k] + a2*dat[k-1] + sin(2*pi*f*k) + rnorm(1)*0.05
}

N = 100
Q=100;
L=rep(10,3)
Pc=0.8
Pm=0.03
max=rep(10,3);min=c(10^-3,10^-3,10^-6)
num_gen=100
lag=2

#test set
k0=53;N_tilda = 100
x_test = matrix(,lag,N_tilda)
for (i in k0:(N_tilda+k0-1)){
  x_test[,i-k0+1] = dat[(i-lag+1):(i)]
}

# output from the commented codes below
load("out.RData")
load("MSE.RData")

#plot step-ahead predictions
step=25
y_test = dat[(k0+step):(k0+N_tilda+step-1)]
plot(y_test,xlab="step k",ylab=paste("y(k+",step,")"))
lines(y_test)
lines(out[step,],col=2);points(out[step,],col=2,pch=2)
legend("topright",c("predictive mean","measurements"),lwd=c(1,1),col=c(2,1),pch=c(2,1))

#perform one to M-step ahead prediction and calculate the MSE
M=30
# MSE=c();out = matrix(,M,N_tilda)
# for (i in 1:3) {
#   res = MGP(dat,N,Q,num_gen,lag,max,min,Pc,Pm,x_test,i)
#   out[i,] = res
#   y_test = dat[(k0+i):(k0+N_tilda+i-1)]
#   MSE[i] = sum((out[i,]-y_test)^2/N_tilda)
#   print(i)
# }

par(mfrow=c(3,3))
for (i in 1:9){
  y_test = dat[(k0+i):(k0+N_tilda+i-1)]
  plot(y_test,xlab="step k",ylab=paste("y(k+",i,")"),ylim=c(-3,3));lines(y_test)
  lines(out[i,],col=2);points(out[i,],col=2,pch=2)
  legend("topright",c("predictive mean","measurements"),lwd=c(1,1),col=c(2,1),pch=c(2,1))
}

par(mfrow=c(1,1))
plot(MSE,xlab="step ahead",ylab="MSE");lines(MSE)

#plot one to M step ahead predictions starting from the (j+k0+1)th entry
j=79
plot((k0+j):(k0+M + j-1),dat[(k0+j):(k0+M + j-1)],xlab="step k",ylab="y(k)",main="M-step ahead prediction");
lines((k0+j):(k0+M + j-1),dat[(k0+j):(k0+M + j-1)])
lines((k0+j):(k0+M + j-1),out[,j][1:30],col=2)
points((k0+j):(k0+M + j-1),out[,j][1:30],col=2,pch=2)
legend("topright",c("predictive mean","measurements"),col=c(2,1),lwd=c(1,1),pch=c(2,1))

