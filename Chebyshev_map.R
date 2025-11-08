library(plot3D)

source("Functions.r")
source("filled.contour3.r")

Chevy <- function(a,n=200){
  x <- rep(NA,n)
  x[1]=0.3
  for (i in 2:n){
    x[i]=cos((a^2)*acos(x[i-1]))
  }
  return(x)
}

Z <- Chevy(0.15,200) 
FCRE(Z,3)

n = 1000
a <- seq(0.1,5,0.05)
q <- seq(0.1,3,0.05)
A <- matrix(NA,length(a),length(q))

for(i in 1:length(a)) for(j in 1:length(q)){
	Z <- Chevy(a[i],n) 
	A[i,j] = FCRE(Z,q[j])
}

image2D(A,a,q,xlab="a",ylab="q")

q=1
n = 1000
a <- seq(0.1,5,0.05)
A <- matrix(NA,length(a),length(a))

for(i in 1:length(a)) for(j in 1:length(a)){
	Z1 <- Chevy(a[i],n) 
	Z2 <- Chevy(a[j],n)
	A[i,j] = FCRI(Z1,Z2,q) #/10
}

image2D(A,a,a,main="q = 1",xlab=expression(a[1]),ylab=expression(a[2]))


q=2
alpha=0.5
n = 1000
a0 = 1
a <- seq(1,5,0.05)
A <- matrix(NA,length(a),length(a))

for(i in 1:length(a)) for(j in 1:length(a)){
	Z1 <- Chevy(a[i],n) 
	Z2 <- Chevy(a[j],n)
	Z3 <- Chevy(a0,n)	
	A[i,j] = JFCRI(Z1,Z2,Z3,q,alpha) #/10
}

image2D(A,a,a,main="q = 2",xlab=expression(a[1]),ylab=expression(a[2]))


### Specific values of a

c = c(0.5,1,2)
n = 1000
q <- seq(0.1,3,0.05)
A <- matrix(NA,length(c),length(q))

for(i in 1:length(c)) for(j in 1:length(q)){
	Z <- Chevy(c[i],n) 
	A[i,j] = FCRE(Z,q[j])
}

dim(A)

plot(q,A[1,],type="l",lwd=2,col="black",ylab="FCRE",ylim=c(0,500))
lines(q,A[2,],col="blue",lwd=2)
lines(q,A[3,],col="red",lwd=2)
legend("top", c("a=0.5","a=1","a=2"),
col=c("black","blue","red"),lwd=2,bty = "n")

