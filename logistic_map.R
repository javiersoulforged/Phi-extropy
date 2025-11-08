library(plot3D)
source("Functions.r")

logistic <- function(c,n=200){
  x <- rep(NA,n)
  x[1]=0.1
  for (i in 2:n){
    x[i] = c * x[i-1] * (1 - x[i-1])
  }
  return(x)
}

Z <- logistic(3,200) 
FCRE(Z,3)

n = 1000
c <- seq(0,4,0.05)
q <- seq(0.1,3,0.05)
A <- matrix(NA,length(c),length(q))

for(i in 1:length(c)) for(j in 1:length(q)){
	Z <- logistic(c[i],n) 
	A[i,j] = FCRE(Z,q[j])
}

image2D(A,c,q,xlab="c",ylab="q")

q=3
n = 1000
c <- seq(0,4,0.05)
A <- matrix(NA,length(c),length(c))

for(i in 1:length(c)) for(j in 1:length(c)){
	Z1 <- logistic(c[i],n) 
	Z2 <- logistic(c[j],n)
	A[i,j] = FCRI(Z1,Z2,q) #/10
}

image2D(A,c,c,main="q = 3",xlab=expression(c[1]),ylab=expression(c[2]))


q=0.5
alpha=0.5
n = 1000
c0 = 4
c <- seq(2,4,0.05)
A <- matrix(NA,length(c),length(c))

for(i in 1:length(c)) for(j in 1:length(c)){
	Z1 <- logistic(c[i],n) 
	Z2 <- logistic(c[j],n)
	Z3 <- logistic(c0,n)	
	A[i,j] = JFCRI(Z1,Z2,Z3,q,alpha) #/10
}

image2D(A,c,c,main="q = 0.5",xlab=expression(c[1]),ylab=expression(c[2]))


### Specific values of c

c = c(3,3.3,3.4)
n = 1000
q <- seq(0.1,3,0.05)
A <- matrix(NA,length(c),length(q))

for(i in 1:length(c)) for(j in 1:length(q)){
	Z <- logistic(c[i],n) 
	A[i,j] = FCRE(Z,q[j])
}

dim(A)

plot(q,A[1,],type="l",lwd=2,col="black",ylab="FCRE",ylim=c(200,500))
lines(q,A[2,],col="blue",lwd=2)
lines(q,A[3,],col="red",lwd=2)
legend("top", c("c=3","c=3.3","c=3.4"),
col=c("black","blue","red"),lwd=2,bty = "n")







